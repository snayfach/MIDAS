#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, shutil
from midas import utility
from midas.merge import merge
from time import time

class GenomicSite:
	def __init__(self, values):
	
		# initialize
		self.ref_id, self.ref_pos, self.ref_allele  = values[0].rsplit('|', 2)
		self.ref_pos = int(self.ref_pos)
		self.id = values[0].rsplit('|', 1)[0]

		# per-sample statistics
		self.sample_counts = [[int(j) for j in i.split(',')] for i in values[1:]] # <list>, per sample allele counts
		self.sample_mafs = []  # <list>, minor allele frequencies
		self.sample_depths = [] # <list>, count of major+minor allele frequencies
		
		# pooled statistics
		self.total_samples = len(self.sample_counts)
		self.pooled_counts = self.compute_pooled_counts() # <list>, count of each allele across all samples
		self.pooled_depth = sum(self.pooled_counts)
		self.major_allele = ''
		self.minor_allele = ''
		self.minor_freq = 0.0
		self.misc_freq = 0.0
					
		# site annotations
		self.site_type = ''
		self.gene_id = ''
		self.amino_acids = ''
		self.ref_codon = ''
		self.codon_pos = ''
		
	def compute_pooled_counts(self):
		""" Compute count of 4 nucleotides across samples """
		pooled_counts = []
		for i in range(4):
			pooled_counts.append(sum([counts[i] for counts in self.sample_counts]))
		return pooled_counts

	def call_alleles(self):
		""" Call major and minor alleles at GenomicSite """
		alleles = list('ATCG')
		pooled_counts = list(self.pooled_counts)

		# major allele
		if max(pooled_counts) == 0: return
		self.major_count = max(pooled_counts)
		self.major_index = pooled_counts.index(self.major_count)
		self.major_allele = alleles[self.major_index]
		alleles.remove(self.major_allele)
		pooled_counts.remove(self.major_count)
	
		# minor allele
		if max(pooled_counts) == 0: return
		self.minor_count = max(pooled_counts)
		self.minor_index = pooled_counts.index(self.minor_count)
		self.minor_allele = alleles[self.minor_index]
		pooled_counts.remove(self.minor_count)
		self.minor_freq = float(self.minor_count)/(self.minor_count+self.major_count)
		
		# tri and quad alleles
		if max(pooled_counts) == 0: return
		self.misc_count = sum(pooled_counts)
		self.misc_freq = self.misc_count/float(self.pooled_depth)

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

	def flag(self, min_maf=0.0, min_prev=0.0, snp_type='any'):
		""" Filter genomic site based on MAF and prevalence """
		if self.minor_freq < min_maf:
			self.flag = (True, 'min_maf')
		elif self.prevalence < min_prev:
			self.flag = (True, 'min_prev')
		elif snp_type == 'bi' and self.misc_freq >= min_maf:
			self.flag = (True, 'snp_type')
		else:
			self.flag = (False, None)

	def annotate(self, genes):
		""" Annotate variant and reference site """
		# genes: list of genes, each gene contains info
		# contig: contig sequence
		# gene_index: current position in list of genes; global variable
		self.amino_acids = {}
				
		while True:
			# 1. fetch next gene
			#    if there are no more genes, snp must be non-coding so break
			if genes[1] < len(genes[0]):
				gene = genes[0][genes[1]]
			else:
				self.site_type = 'NC'
				return
			# 2. if snp is upstream of next gene, snp must be non-coding so break
			if (self.ref_id < gene['scaffold_id'] or
			   (self.ref_id == gene['scaffold_id'] and self.ref_pos < gene['start'])):
				self.site_type = 'NC'
				return
			# 3. if snp is downstream of next gene, pop gene, check (1) and (2) again
			if (self.ref_id > gene['scaffold_id'] or
			   (self.ref_id == gene['scaffold_id'] and self.ref_pos > gene['end'])):
				genes[1] += 1
				continue
			# 4. otherwise, snp must be in gene
			#    annotate site (1D-4D)
			else:
				self.gene_id = gene['gene_id']
				self.ref_codon, self.codon_pos = self.fetch_ref_codon(gene)
				if not all([_ in ['A','T','C','G'] for _ in self.ref_codon]): # check for invalid bases in codon
					self.site_type = 'NA'
				else:
					for allele in ['A','T','C','G']: # + strand
						codon = utility.index_replace(self.ref_codon, allele, self.codon_pos, gene['strand']) # +/- strand
						self.amino_acids[allele] = utility.translate(codon)
					unique_aa = set(self.amino_acids.values())
					degeneracy = 4 - len(unique_aa) + 1
					self.site_type = '%sD' % degeneracy
					# AA's identical: degeneracy = 4 - 1 + 1 = 4
					# AA's all different, degeneracy = 4 - 4 + 1 = 1
				return

	def fetch_ref_codon(self, gene):
		""" Fetch codon within gene for given site """
		# position of site in gene
		gene_pos = self.ref_pos-gene['start'] if gene['strand']=='+' else gene['end']-self.ref_pos
		# position of site in codon
		codon_pos = gene_pos % 3
		# gene sequence (oriented start to stop)
		ref_codon = gene['seq'][gene_pos-codon_pos:gene_pos-codon_pos+3]
		return ref_codon, codon_pos
						
	def write(self, files):
		""" Store data for GenomicSite in Species"""
		# snps_info
		atcg_counts = ','.join([str(_) for _ in self.pooled_counts])
		atcg_aas = ','.join([self.amino_acids[_] for _ in list('ATCG')]) if self.amino_acids != {} else ''
		info = [self.id,
				self.ref_id,
				str(self.ref_pos),
				self.ref_allele,
				self.major_allele,
				self.minor_allele,
				'{0:.3g}'.format(self.minor_freq),
				str(self.count_samples),
				atcg_counts,
				self.site_type,
				atcg_aas,
				self.gene_id]
		info = '\t'.join(info)+'\n'
		files['info'].write(info)
		# snps_freq
		freq = self.id + '\t' + '\t'.join(['{0:.3g}'.format(freq) for freq in self.sample_mafs])+'\n'
		files['freq'].write(freq)
		# snps_depth
		depth = self.id + '\t' + '\t'.join([str(depth) for depth in self.sample_depths])+'\n'
		files['depth'].write(depth)

def read_run_midas_snps(species_id, samples):
	""" Open SNP files for species across samples """
	infiles = []
	for sample in samples:
		path = '%s/snps/output/%s.snps.gz' % (sample.dir, species_id)
		file = utility.iopen(path)
		next(file)
		infiles.append(file)
	return infiles

def write_count_matrix(outdir, sample_ids, thread=None):
	""" Open matrix for writing and write header """
	if thread is None:
		file = open('%s/atcg_counts.txt' % outdir, 'w')
	else:
		file = open('%s/atcg_counts.%s.txt' % (outdir, thread), 'w')
	file.write('\t'.join(['site_id']+sample_ids)+'\n')
	return file

def build_temp_count_matrix(tempdir, species_id, samples, thread, max_sites):
	""" Build SNP matrices using a subset of total samples """
	sample_ids = [s.id for s in samples]
	midas_files = read_run_midas_snps(species_id, samples)
	matrix_file = write_count_matrix(tempdir, sample_ids, thread)
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
			allele_counts = [r[-1] for r in records]
			matrix_file.write(site_id+'\t'+'\t'.join(allele_counts)+'\n')

def parallel(function, argument_list, threads):
	""" Based on: https://gist.github.com/admackin/003dd646e5fadee8b8d6 """
	import multiprocessing as mp
	import signal
	import time
	
	def init_worker():
		signal.signal(signal.SIGINT, signal.SIG_IGN)

	pool = mp.Pool(threads, init_worker)
	
	try:
		results = []
		for arguments in argument_list:
			p = pool.apply_async(function, args=arguments)
			results.append(p)
		pool.close()
		while True:
			if all(r.ready() for r in results):
				return
			time.sleep(1)

	except KeyboardInterrupt:
		pool.terminate()
		pool.join()
		sys.exit("\nKeyboardInterrupt")

def parallel_build_temp_count_matrixes(species, args):
	""" Split up samples into batches, merge each batch, merge together batches """
	samples_list = utility.batch_samples(species.samples, threads=args['threads'])	
	argument_list = []
	for thread, sample_ids in zip(range(args['threads']), samples_list):
		arguments=(species.tempdir, species.id, sample_ids, thread, args['max_sites'])
		argument_list.append(arguments)
	parallel(build_temp_count_matrix, argument_list, args['threads'])

def	read_count_matrixes(species, args):
	""" Open matrices for reading and skip headers """
	files = []
	for j in range(args['threads']):
		path = '%s/atcg_counts.%s.txt' % (species.tempdir, j)
		files.append(open(path))
		next(files[j])
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
	info_fields = ['site_id', 'ref_id', 'ref_pos',
				   'ref_allele', 'major_allele', 'minor_allele', 'minor_freq',
				   'count_samples', 'count_atcg', 'site_type',
				   'amino_acid_atcg', 'gene_id']
	files['info'].write('\t'.join(info_fields)+'\n')
	return files

def build_sharded_tables(species, args, thread):
	""" Build merged output files for species using every <thread>th site"""
	infiles = read_count_matrixes(species, args)
	outfiles = write_merge_midas(species, args, thread)
	genes = utility.read_genes(species.id, args['db'])
	line_num = -1
	while True:
		# fetch allele counts for next site
		try:
			lines = [next(file) for file in infiles]
		except StopIteration:
			break
		# control which lines are processed
		line_num += 1
		if (line_num % args['threads'] == thread) is False:
			continue
		# initialize GenomicSite
		values = sum([line.rstrip().split()[1:] for line in lines], [lines[0].split()[0]])
		site = GenomicSite(values)
		site.call_alleles()
		site.compute_per_sample_mafs()
		site.compute_prevalence(species.sample_depth, args['site_depth'], args['site_ratio'])
		site.flag(args['snp_maf'], args['site_prev'], args['snp_type'])
		# decide what to do with site
		if site.flag[0] is True:
			continue
		else:
			site.annotate(genes)
			site.write(outfiles)
	# finish up
	for file in infiles: file.close()
	for file in outfiles.values(): file.close()

def parallel_build_sharded_tables(species, args):
	""" """
	argument_list = []
	for thread in range(args['threads']):
		arguments=(species, args, thread)
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
	# write headers
	for thread in range(args['threads']):
		for ftype in ['freq', 'depth', 'info']:
			next(infiles[thread][ftype])
	# write remainder of lines
	while True:
		for thread in range(args['threads']):
			try:
				outfiles['freq'].write(next(infiles[thread]['freq']))
				outfiles['depth'].write(next(infiles[thread]['depth']))
				outfiles['info'].write(next(infiles[thread]['info']))
			except StopIteration:
				return

def run_pipeline(args):
	
	print("Identifying species and samples")
	species_list = merge.select_species(args, dtype='snps')
	for species in species_list:
		print("  %s" % species.id)
		print("    genome name: %s" % species.genome_info['genome_name'])
		print("    genome length: %s" % species.genome_info['length'])
		print("    count contigs: %s" % max(1, int(species.genome_info['contigs'])))
		print("    count samples: %s" % len(species.samples))
	
	print("\nMerging snps")
	for species in species_list:
	
		print("  %s" % species.id)
		species.tempdir = '%s/%s/temp' % (args['outdir'], species.id)
		if not os.path.isdir(species.tempdir): os.mkdir(species.tempdir)

		print("    merging count data")
		parallel_build_temp_count_matrixes(species, args)

		print("    calling SNPs")
		parallel_build_sharded_tables(species, args)

		print("    writing output files")
		merge_sharded_tables(species, args)

		print("    finishing")
		merge.write_snps_readme(args, species)
		species.write_sample_info(dtype='snps', outdir=args['outdir'])
		shutil.rmtree(species.tempdir)
		print "    done!"


