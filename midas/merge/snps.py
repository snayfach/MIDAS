#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, shutil
from midas import utility
from midas.merge import merge
from time import time
from operator import itemgetter
import multiprocessing
from smelter.iggdb import IGGdb
from smelter.utilities import tserr
import traceback

class GenomicSite:
	def __init__(self, id, values):
	
		# initialize
		self.id = str(id) #values[0].rsplit('|', 1)[0]
		self.ref_id, self.ref_pos, self.ref_allele  = values[0].rsplit('|', 2)
		self.ref_pos = int(self.ref_pos)
		
		# per-sample statistics
		self.sample_counts = [[int(j) for j in i.split(',')] for i in values[1:]] # <list>, per sample allele counts
		self.sample_mafs = []  # <list>, minor allele frequencies
		self.sample_depths = [] # <list>, count of major+minor allele frequencies
		
		# pooled statistics
		self.total_samples = len(self.sample_counts)
		self.pooled_counts = self.compute_pooled_counts() # <list>, count of each allele across all samples
		self.pooled_depth = sum(self.pooled_counts)
		self.major_allele = None
		self.minor_allele = None
		self.snp_type = None # mono, bi, tri, quad
					
		# site annotations
		self.locus_type = None # CDS, tRNA, rRNA, IGR
		self.site_type = None # 1D, 2D, 3D, 4D
		self.gene_id = None
		self.amino_acids = None
		self.ref_codon = None
		self.codon_pos = None
		
	def compute_pooled_counts(self):
		""" Compute count of 4 nucleotides across samples """
		pooled_counts = []
		for i in range(4):
			pooled_counts.append(sum([counts[i] for counts in self.sample_counts]))
		return pooled_counts

	def call_alleles(self, snp_freq):
		""" Call major and minor alleles at GenomicSite """
		
		# check at least 1 mapped read
		if self.pooled_depth == 0: return
		
		# get sorted counts
		alleles = list('ACGT')
		freqs = [float(count)/self.pooled_depth for count in self.pooled_counts]
		allele_freqs = sorted(zip(alleles, freqs), key=itemgetter(1), reverse=True)

		# major allele
		if allele_freqs[0][1] > 0:
			self.major_allele = allele_freqs[0][0]
			self.major_index = alleles.index(self.major_allele)
		
		# minor allele
		if allele_freqs[1][1] > 0:
			self.minor_allele = allele_freqs[1][0]
			self.minor_index = alleles.index(self.minor_allele)
		
		# classify SNP
		snp_types = ['quad', 'tri', 'bi', 'mono']
		for snp_type, allele_freq in zip(snp_types, allele_freqs[::-1]):
			allele, freq = allele_freq
			if freq >= snp_freq:
				self.snp_type = snp_type
				break

	def compute_per_sample_mafs(self):
		""" Compute per-sample depth (major + minor allele) and minor allele freq at GenomicSite """
		if self.major_allele is None:
			self.sample_mafs = [0.0 for counts in self.sample_counts]
			self.sample_depths = [0 for counts in self.sample_counts]
		elif self.minor_allele is None:
			self.sample_mafs = [0.0 for counts in self.sample_counts]
			self.sample_depths = [counts[self.major_index] for counts in self.sample_counts]
		else:
			for counts in self.sample_counts:
				sample_depth = counts[self.major_index] + counts[self.minor_index]
				sample_maf = float(counts[self.minor_index])/sample_depth if sample_depth > 0 else 0.0
				self.sample_mafs.append(sample_maf)
				self.sample_depths.append(sample_depth)

	def compute_prevalence(self, mean_depths, min_depth, max_ratio):
		""" Compute the fraction of samples where site passes all filters """
		pass_qc = []
		for mean_depth, site_depth in zip(mean_depths, self.sample_depths):
			if site_depth < min_depth:
				pass_qc.append(0)
			elif site_depth/mean_depth > max_ratio:
				pass_qc.append(0)
			else:
				pass_qc.append(1)
		self.count_samples = sum(pass_qc)
		self.prevalence = sum(pass_qc)/float(len(pass_qc))

	def flag(self, min_prev, snp_types):
		""" Filter genomic site based on MAF and prevalence """
		if self.prevalence < min_prev:
			self.flag = (True, 'min_prev')
		elif ('any' not in snp_types and 
				self.snp_type not in snp_types):
			self.flag = (True, 'snp_type')
		else:
			self.flag = (False, None)

	def annotate(self, genes):
		""" Annotate variant and reference site """
		# genes = {'list': list of genes, 'index':index position in list}
		# each element in 'list' is sorted by scaffold_id (asc), start (asc), end (desc)
		# each element is a dictionary with gene info
		while True:
			# 1. fetch next gene
			#    if there are no more genes, snp must be intergenic so break
			if genes['index'] < len(genes['list']):
				gene = genes['list'][genes['index']]
			else:
				self.locus_type = 'IGR'
				return
			# 2. if snp is upstream of next gene, snp must be intergenic so break
			if (self.ref_id < gene['scaffold_id'] or
			   (self.ref_id == gene['scaffold_id'] and self.ref_pos < gene['start'])):
				self.locus_type = 'IGR'
				return
			# 3. if snp is downstream of next gene, pop gene, check (1) and (2) again
			if (self.ref_id > gene['scaffold_id'] or
			   (self.ref_id == gene['scaffold_id'] and self.ref_pos > gene['end'])):
				genes['index'] += 1
				continue
			# 4. snp in coding gene: annotate (1D-4D)
			elif gene['gene_type'] == 'CDS':
				self.locus_type = gene['gene_type']
				self.gene_id = gene['gene_id']
				if len(gene['seq']) % 3 != 0: # gene must by divisible by 3 to id codons
					return
				self.ref_codon, self.codon_pos = self.fetch_ref_codon(gene)
				if not all([_ in ['A','T','C','G'] for _ in self.ref_codon]): # codon can't contain weird characters
					return
				self.amino_acids = []
				for allele in ['A','C','G','T']: # + strand
					codon = utility.index_replace(self.ref_codon, allele, self.codon_pos, gene['strand']) # +/- strand
					amino_acid = utility.translate(codon)
					self.amino_acids.append(amino_acid)
				unique_aa = set(self.amino_acids)
				degeneracy = 4 - len(unique_aa) + 1
				self.site_type = '%sD' % degeneracy
				self.amino_acids = ','.join(self.amino_acids)
				# AA's identical: degeneracy = 4 - 1 + 1 = 4
				# AA's all different, degeneracy = 4 - 4 + 1 = 1
				return
			# 5. snp in non-coding gene
			else:
				self.locus_type = gene['gene_type']
				self.gene_id = gene['gene_id']
				return

	def fetch_ref_codon(self, gene):
		""" Fetch codon within gene for given site """
		# position of site in gene
		gene_pos = self.ref_pos - gene['start'] if gene['strand'] == '+' else gene['end'] - self.ref_pos
		# position of site in codon
		codon_pos = gene_pos % 3
		# gene sequence (oriented start to stop)
		ref_codon = gene['seq'][gene_pos-codon_pos:gene_pos-codon_pos+3]
		return ref_codon, codon_pos
						
	def write(self, files):
		""" Store data for GenomicSite in Species"""
		# snps_info
		count_a, count_c, count_g, count_t = [str(_) for _ in self.pooled_counts]
		info = [self.id,
				self.ref_id,
				str(self.ref_pos),
				self.ref_allele,
				self.major_allele,
				self.minor_allele,
				str(self.count_samples),
				count_a, count_c, count_g, count_t,
				self.locus_type,
				self.gene_id,
				self.snp_type,
				self.site_type,
				self.amino_acids,
				]
		info = '\t'.join([replace_none(_) for _ in info])+'\n'
		files['info'].write(info)
		# snps_freq
		freq = self.id + '\t' + '\t'.join(['{0:.3g}'.format(freq) for freq in self.sample_mafs])+'\n'
		files['freq'].write(freq)
		# snps_depth
		depth = self.id + '\t' + '\t'.join([str(depth) for depth in self.sample_depths])+'\n'
		files['depth'].write(depth)


def parallel(function, argument_list, threads):
	mp = multiprocessing.Pool(threads)
	return mp.starmap(function, argument_list)


def replace_none(input_string, replace_string="NA"):
	return input_string if input_string is not None else replace_string

def read_run_midas_snps(species_id, samples):
	""" Open SNP files for species across samples """
	infiles = []
	for sample in samples:
		path = '%s/snps/output/%s.snps.gz' % (sample.dir, species_id)
		file = utility.iopen(path)
		next(file)
		infiles.append(file)
	return infiles

def build_temp_count_matrix(tempdir, species_id, samples, split_num, max_sites):
	""" Build SNP matrices using a subset of total samples """
	sample_ids = [s.id for s in samples]
	midas_files = read_run_midas_snps(species_id, samples)
	matrix_file = open('%s/acgt_counts.%s.txt' % (tempdir, split_num), 'w')
	matrix_file.write('\t'.join(['site_id']+sample_ids)+'\n')
	nsites = 0
	while True:
		records = []
		for file in midas_files:
			try:
				record = next(file).rstrip('\n').split('\t')
				records.append(record)
			except StopIteration:
				matrix_file.close()
				for file in midas_files: file.close()
				return
		if nsites >= max_sites:
			matrix_file.close()
			for file in midas_files: file.close()
			return
		else:
			nsites += 1
			site_id = '|'.join(records[0][0:3])
			allele_counts = [','.join(r[-4:]) for r in records]
			matrix_file.write(site_id+'\t'+'\t'.join(allele_counts)+'\n')

def parallel_build_temp_count_matrixes(species, args):
	""" Split up samples into batches, merge each batch, merge together batches """
	argument_list = []
	for split_num, sample_ids in enumerate(species.sample_lists):
		arguments=(species.tempdir, species.id, sample_ids, split_num, args['max_sites'])
		argument_list.append(arguments)
	parallel(build_temp_count_matrix, argument_list, args['threads'])

def	read_count_matrixes(species, args):
	""" Open matrices for reading and skip headers """
	files = []
	for split_num in range(species.num_splits):
		path = '%s/acgt_counts.%s.txt' % (species.tempdir, split_num)
		files.append(open(path))
		next(files[split_num])
	return files

def write_merge_midas(species, args, thread=None):
	""" Open output files for species """
	files = {}
	# open files
	for ftype in ['info', 'freq', 'depth']:
		if thread is None:
			path = '%s/%s/snps_%s.txt' % (args['outdir'], species.id, ftype)
		else:
			path = '%s/%s/temp/snps_%s.%s.txt' % (args['outdir'], species.id, ftype, thread)
		files[ftype] = open(path, 'w')
	# write headers
	for ftype in ['freq', 'depth']:
		record = ['site_id']+[s.id for s in species.samples]
		files[ftype].write('\t'.join(record)+'\n')
	info_fields = ['site_id', 
				   'ref_id', 
				   'ref_pos',
				   'ref_allele', 
				   'major_allele', 
				   'minor_allele',
				   'count_samples', 
				   'count_a', 
				   'count_c',
				   'count_g',
				   'count_t',
				   'locus_type',
				   'gene_id',
				   'snp_type',
				   'site_type',
				   'amino_acids'
				   ]
	files['info'].write('\t'.join(info_fields)+'\n')
	return files

def build_sharded_tables(species, args, thread, line_from, line_to, errcnt=multiprocessing.Value('i', 0)):
	""" Build merged output files for species using every <thread>th site"""
	infiles = read_count_matrixes(species, args)
	outfiles = write_merge_midas(species, args, thread)
	try:
		genes = utility.read_genes(species.id, args['iggdb'])
	except:
		genes = None
		with errcnt.get_lock():
			errcnt.value += 1
			if errcnt.value == 1:
				traceback.print_exc()
				tserr("Apologies - gene annotations disabled in this run.")
	
	line_num = -1
	while True:
				
		# fetch allele counts for next site
		try:
			lines = [next(file) for file in infiles]
		except StopIteration:
			break
		
		# control which lines are processed
		line_num += 1
		if line_num < line_from:
			continue
		elif line_num >= line_to:
			break
		
		# initialize GenomicSite
		values = sum([line.rstrip().split()[1:] for line in lines], [lines[0].split()[0]])
		site_id = line_num+1
		site = GenomicSite(site_id, values)
		site.call_alleles(args['allele_freq'])
		site.compute_per_sample_mafs()
		site.compute_prevalence(species.sample_depth, args['site_depth'], args['site_ratio'])
		site.flag(args['site_prev'], args['snp_type'])
		
		# decide what to do with site
		if site.flag[0] is True:
			continue
		if genes:
			site.annotate(genes)
		site.write(outfiles)
	
	# finish up
	for file in infiles: file.close()
	for file in outfiles.values(): file.close()

def parallel_build_sharded_tables(species, args):
	""" """
	# get number of total lines to process
	infiles = read_count_matrixes(species, args)
	num_lines = 0
	for line in infiles[0]:
		num_lines += 1
	for file in infiles:
		file.close()

	lines_per = max(1, num_lines/args['threads'])
	line_ranges = [[thread * lines_per, thread * lines_per + lines_per] for thread in range(args['threads'])]
	line_ranges[-1][-1] = num_lines

	argument_list = []
	for thread, line_range in enumerate(line_ranges):
		line_from, line_to = line_range
		arguments=(species, args, thread, line_from, line_to)
		argument_list.append(arguments)
	
	parallel(build_sharded_tables, argument_list, args['threads'])

def merge_sharded_tables(species, args):
	""" Merge N sets of sharded tables, where N is the number of threads"""
	# open output tables
	outfiles = write_merge_midas(species, args)
	# open input sharded tables
	infiles = {}
	for thread in range(args['threads']):
		infiles[thread] = {}
		for ftype in ['freq', 'depth', 'info']:
			path = '%s/%s/temp/snps_%s.%s.txt' % (args['outdir'], species.id, ftype, thread)
			infiles[thread][ftype] = open(path)
	# skip headers
	for thread in range(args['threads']):
		for ftype in ['freq', 'depth', 'info']:
			next(infiles[thread][ftype])
	# write remainder of lines
	for thread in range(args['threads']):
		for ftype in ['freq', 'depth', 'info']:
			for line in infiles[thread][ftype]:
				outfiles[ftype].write(line)

def write_snps_readme(args, sp):
	outfile = open('%s/%s/readme.txt' % (args['outdir'], sp.id), 'w')
	outfile.write("""
Description of output files and file formats from 'merge_midas.py snps'

Output files
############
snps_freq.txt
  frequency of minor allele per genomic site and per sample
  a value of 1.0 indicates that all reads matched the minor allele for site-sample
  the major (most common) and minor allele (2nd most common) are determined from pooled reads across ALL samples
  see: snps_info.txt for details on the major, minor, and reference alleles
snps_depth.txt
  number of reads mapped to genomic site per sample
  only accounts for reads matching either major or minor allele
snps_info.txt  
  metadata for genomic site
  see below for more information
snps_summary.txt
  alignment summary statistics per sample
  see below for more information
snps_log.txt
  log file containing parameters used

Output formats
############
snps_freq.txt and snps_depth.txt
  tab-delimited matrix files
  field names are sample ids
  row names are genome site ids
  see: snps_info.txt for details on each genomic site
snps_summary.txt
  sample_id: sample identifier
  genome_length: number of base pairs in representative genome
  covered_bases: number of reference sites with at least 1 mapped read
  fraction_covered: proportion of reference sites with at least 1 mapped read
  mean_coverage: average read-depth across reference sites with at least 1 mapped read
  aligned_reads: number of reads that aligned BEFORE quality filtering
  mapped_reads: number of reads that aligned AFTER quality filtering
snps_info.txt
  site_id: incrementing integer field
  ref_id: identifier of scaffold in representative genome
  ref_pos: position of site on ref_id
  ref_allele: allele in reference genome
  major_allele: most common allele in metagenomes
  minor_allele: second most common allele in metagenomes
  count_samples: number of metagenomes where site_id was found
  count_a: count of A allele in pooled metagenomes
  count_c: count of C allele in pooled metagenomes
  count_g: count of G allele in pooled metagenomes
  count_t: count of T allele in pooled metagenomes
  locus_type: CDS (site in coding gene), RNA (site in non-coding gene), IGR (site in intergenic region)
  gene_id: gene identified if locus_type is CDS, or RNA
  snp_type: indicates the number of alleles observed at site (mono,bi,tri,quad); observed allele are determined by --snp_maf flag  
  site_type: indicates degeneracy: 1D, 2D, 3D, 4D
  amino_acids: amino acids encoded by 4 possible alleles

Additional information for species can be found in the reference database:
 %s/rep_genomes/%s
""" % (args['db'], sp.id) )
	outfile.close()


def per_species_work(index):
	global species_list
	global global_args
	args = global_args
	species = species_list[index]
	print("  %s" % species.id)
	species.tempdir = '%s/%s/temp' % (args['outdir'], species.id)
	if not os.path.isdir(species.tempdir): os.mkdir(species.tempdir)
	species.sample_lists = utility.batch_samples(species.samples, threads=args['threads'])
	species.num_splits = len(species.sample_lists)
	
	print("    merging count data")
	parallel_build_temp_count_matrixes(species, args)

	print("    calling SNPs")
	parallel_build_sharded_tables(species, args)

	# this is the slow step, and it is single-threaded
	# that is why we have to run this function in multiprocessing,
	# despite the fact the other steps above are "parallel"
	print("    writing output files")
	merge_sharded_tables(species, args)

	print("    finishing")
	write_snps_readme(args, species)
	species.write_sample_info(dtype='snps', outdir=args['outdir'])
	shutil.rmtree(species.tempdir)


def psw_safe(index, sem):
	try:
		per_species_work(index)
	finally:
		sem.release()

	
def run_pipeline(args):
	
	print("Identifying species and samples")
	if 'db' in args:
		args['iggdb'] = IGGdb(f"{args['db']}/metadata/species_info.tsv")
	global species_list
	species_list = merge.select_species(args, dtype='snps')
	for species in species_list:
		print("  %s" % species.id)
		if 'genome_name' in species.genome_info: 
			print("    genome name: %s" % species.genome_info['genome_name'])
		if 'length' in species.genome_info:
			print("    genome length: %s" % species.genome_info['length'])
		if 'contigs' in species.genome_info:
			print("    count contigs: %s" % max(1, int(species.genome_info['contigs'])))
		print("    count samples: %s" % len(species.samples))
	
	print("\nMerging snps")
	
	global global_args
	global_args = args
	sem = multiprocessing.Semaphore(int(args['threads']))
	procs = []
	for index in range(0, len(species_list)):
		sem.acquire()
		procs.append(multiprocessing.Process(target=psw_safe, args=[index, sem]))
		procs[-1].start()
	for p in procs:
		p.join()

