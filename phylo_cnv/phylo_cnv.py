#!/usr/bin/python

# PhyloCNV - estimating the abundance, gene-content, and phylogeny of microbes from metagenomes
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.1'

# Libraries
# ---------
import sys
import os
import numpy as np
import argparse
import pysam
import gzip
import time
import subprocess
import operator
import Bio.SeqIO
import phylo_species
import resource
from collections import defaultdict
from math import ceil
from multiprocessing import Process

# Functions
# ---------

def parallel_process(function, args_list, threads, verbose):
	""" Run function using multiple threads """
	processes = []
	for pargs in args_list: # run function for each set of args in args_list
		p = Process(target=function, kwargs=pargs)
		processes.append(p)
		p.start()
		# control number of active processes
		while len(processes) >= threads:
			indexes = []
			# keep processes that are still alive
			for index, process in enumerate(processes):
				if process.is_alive(): indexes.append(index)
			processes = [processes[i] for i in indexes]
	# wait until there are no active processes
	while len(processes) > 0:
		indexes = []
		for index, process in enumerate(processes):
			if process.is_alive(): indexes.append(index)
		processes = [processes[i] for i in indexes]

def check_arguments(args):
	""" Check validity of command line arguments """
	
	# Pipeline options
	if not any([args['all'], args['species_profile'],
				args['pangenome_build_db'], args['pangenome_align'], args['pangenome_cov'],
				args['snps_build_db'], args['snps_align'], args['snps_call']]):
		sys.exit('Specify one or more pipeline option(s): --all, --profile, --align, --map, --cov, --extract, --remap, --snps')
	if args['all']:
		args['species_profile'] = True
		args['pangenome_build_db'] = True
		args['pangenome_align'] = True
		args['pangenome_cov'] = True
		args['snps_build_db'] = True
		args['snps_align'] = True
		args['snps_call'] = True
	if args['tax_mask'] and not args['tax_map']:
		sys.exit('Specify file mapping read ids in FASTQ file to genome ids in reference database')

	# Input options
	if not args['m1'] and (args['profile'] or args['align']):
		sys.exit('Specify input FASTQ file(s) with -1 -2 or -U')
	if args['m1'] and not os.path.isfile(args['m1']):
		sys.exit('Input file specified with -1 does not exist')
	if args['m2'] and not os.path.isfile(args['m2']):
		sys.exit('Input file specified with -2 does not exist')
	if args['db_dir'] and not os.path.isdir(args['db_dir']):
		sys.exit('Input directory specified with --db-dir does not exist')

	# Output options
	if not args['out']:
		sys.exit('Specify output directory with -o')

def add_binaries(args):
	""" Add paths to external binaries """
	main_dir = os.path.dirname(os.path.abspath(__file__))
	args['bowtie2-build'] = '/'.join([main_dir, 'bin', 'bowtie2-build'])
	args['bowtie2'] = '/'.join([main_dir, 'bin', 'bowtie2'])
	args['samtools'] = '/'.join([main_dir, 'bin', 'samtools'])
	args['bedcov'] = '/'.join([main_dir, 'bin', 'coverageBed'])

def print_copyright():
	# print out copyright information
	print ("-------------------------------------------------------------------------")
	print ("PhyloCNV - estimating the abundance, gene-content, and phylogeny of microbes from metagenomes")
	print ("version %s; github.com/snayfach/PhyloCNV" % __version__)
	print ("Copyright (C) 2015 Stephen Nayfach")
	print ("Freely distributed under the GNU General Public License (GPLv3)")
	print ("-------------------------------------------------------------------------")

def read_phylo_species(inpath):
	""" Parse output from PhyloSpecies """
	if not os.path.isfile(inpath):
		sys.exit("Could not locate species profile: %s\nTry rerunning with --profile" % inpath)
	dict = {}
	fields = [
		('cluster_id', str), ('reads', float), ('bp', float), ('rpkg', float),
		('cov', float), ('prop_cov', float), ('rel_abun', float)]
	infile = open(inpath)
	next(infile)
	for line in infile:
		values = line.rstrip().split()
		dict[values[0]] = {}
		for field, value in zip(fields[1:], values[1:]):
			dict[values[0]][field[0]] = field[1](value)
	return dict

def select_genome_clusters(cluster_abundance, args):
	""" Select genome clusters to map to """
	my_clusters = {}
	# prune all genome clusters that are missing from database
	# this can happen when using an environment specific database
	for cluster_id in cluster_abundance.copy():
		if not os.path.isdir('/'.join([args['db_dir'], cluster_id])):
			del cluster_abundance[cluster_id]
	# user specified a single genome-cluster
	if args['gc_id']:
		cluster_id = args['gc_id']
		if cluster_id not in cluster_abundance:
			sys.exit("Error: specified genome-cluster id %s not found" % cluster_id)
		else:
			abundance = cluster_abundance[args['gc_id']]['rel_abun']
			my_clusters[args['gc_id']] = abundance
	# user specified a list of genome-clusters
	elif args['gc_list']:
		for cluster_id in args['gc_list'].split(','):
			if cluster_id not in cluster_abundance:
				sys.exit("Error: specified genome-cluster id %s not found" % cluster_id)
			else:
				abundance = cluster_abundance[cluster_id]['rel_abun']
				my_clusters[cluster_id] = coverage
	# user specifed a coverage threshold
	elif args['gc_cov']:
		for cluster_id, values in cluster_abundance.items():
			if values['cell_count'] >= args['gc_cov']:
				my_clusters[cluster_id] = values['cov']
	# user specifed a relative-abundance threshold
	elif args['gc_rbun']:
		for cluster_id, values in cluster_abundance.items():
			if values['prop_mapped'] >= args['gc_rbun']:
				my_clusters[cluster_id] = values['rel_abun']
	# user specifed a relative-abundance threshold
	elif args['gc_topn']:
		cluster_abundance = [(i,d['rel_abun']) for i,d in cluster_abundance.items()]
		sorted_abundance = sorted(cluster_abundance, key=operator.itemgetter(1), reverse=True)
		for cluster_id, coverage in sorted_abundance[0:args['gc_topn']]:
			my_clusters[cluster_id] = coverage
	return my_clusters


def pangenome_align(args, tax_mask):
	""" Use Bowtie2 to map reads to all specified genome clusters """
	# Build command
	command = '%s --no-unal ' % args['bowtie2']
	#   index
	command += '-x %s ' % '/'.join([args['out'], 'db', 'pangenomes'])
	#   specify reads
	if args['reads_align']: command += '-u %s ' % args['reads_align']
	#   speed/sensitivity
	command += '--%s ' % args['align_speed']
	#   threads
	command += '--threads %s ' % args['threads']
	#	report up to 20 hits/read if masking hits
	if tax_mask: command += '-k 20 '
	#   input file
	if (args['m1'] and args['m2']): command += '-1 %s -2 %s ' % (args['m1'], args['m2'])
	else: command += '-U %s' % args['m1']
	#   output
	bampath = '/'.join([args['out'], 'pangenome.bam'])
	command += '| %s view -b - > %s' % (args['samtools'], bampath)
	# Run command
	if args['verbose']: print("    running: %s") % command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	#sys.stderr.write(err) # write to stderr: bowtie2 output

def genome_align(args):
	""" Use Bowtie2 to map reads to representative genomes from each genome cluster
	"""
	# Build command
	#	bowtie2
	command = '%s --no-unal ' % args['bowtie2']
	#   speed/sensitivity
	command += '--%s ' % args['align_speed']
	#   threads
	command += '--threads %s ' % args['threads']
	#   bt2 index
	command += '-x %s ' % '/'.join([args['out'], 'db', 'genomes'])
	#   input file
	if (args['m1'] and args['m2']): command += '-1 %s -2 %s ' % (args['m1'], args['m2'])
	else: command += '-U %s' % args['m1']
	#   convert to bam
	command += '| %s view -b - ' % args['samtools']
	#   sort bam
	command += '| %s sort -f - %s ' % (args['samtools'], os.path.join(args['out'], 'genomes.bam'))
	# Run command
	if args['verbose']: print("    running: %s") % command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()

def pileup(args):
	""" Use Samtools to create pileup, filter low quality bases, and write results to VCF file """
	# Build command
	#   mpileup
	command = '%s mpileup -uv -A -d 10000 --skip-indels -B ' % args['samtools']
	#   quality filtering
	command += '-q %s -Q %s ' % (args['snps_mapq'], args['snps_baseq'])
	#   reference fna file
	command += '-f %s ' % '/'.join([args['out'], 'db/genomes.fa'])
	#   input bam file
	command += '%s ' % '/'.join([args['out'], 'genomes.bam'])
	#   output vcf file
	command += '> %s ' % '/'.join([args['out'], 'genomes.vcf'])
	# Run command
	if args['verbose']: print("    running: %s") % command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
		
def format_vcf(args, genome_clusters):
	""" Format vcf output for easy parsing """
	# map scaffold to cluster_id
	ref_to_cluster = {}
	for line in open('/'.join([args['out'], 'db/genomes.map'])):
		ref_id, cluster_id = line.rstrip().split()
		ref_to_cluster[ref_id] = cluster_id
	# open outfiles for each cluster_id
	outdir = '/'.join([args['out'], 'snps'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	outfiles = {}
	for cluster_id in genome_clusters:
		outfiles[cluster_id] = gzip.open('/'.join([outdir, '%s.snps.gz' % cluster_id]), 'w')
		outfiles[cluster_id].write('\t'.join(['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
											  'count_alleles', 'count_ref', 'count_alt', 'depth', 'ref_freq'])+'\n')
	# parse vcf into snp files for each cluster_id
	inpath = '/'.join([args['out'], 'genomes.vcf'])
	for r in parse_vcf(inpath):
		rec = [r['ref_id'], r['ref_pos'], r['ref_allele'], r['alt_allele'], r['cons_allele'],
			   r['count_alleles'], r['count_ref'], r['count_alt'], r['depth'], r['ref_freq']]
		outfile = outfiles[ref_to_cluster[r['ref_id']]]
		outfile.write('\t'.join([str(x) for x in rec])+'\n')

def parse_vcf(inpath):
	""" Yields formatted records from VCF output """
	infile = open(inpath)
	for line in infile:
		# skip header and split line
		if line[0] == '#': continue
		r = line.rstrip().split()
		# get alt alleles
		alt_alleles = r[4].split(',')
		if '<X>' in alt_alleles: alt_alleles.remove('<X>')
		count_alleles = 1 + len(alt_alleles)
		# get allele counts
		info = dict([(_.split('=')) for _ in r[7].split(';')])
		counts = [int(_) for _ in info['I16'].split(',')[0:4]]
		# get consensus allele
		# *note: occassionally there are counts for alternate alleles, but no listed alternate alleles
		if sum(counts) == 0:
			cons_allele = 'NA'
		elif sum(counts[0:2]) >= sum(counts[2:4]):
			cons_allele = r[3]
		elif len(alt_alleles) == 0:
			cons_allele = 'NA'
		else:
			cons_allele = alt_alleles[0]
		# yield formatted record
		yield {'ref_id':r[0],
			   'ref_pos':r[1],
			   'ref_allele':r[3],
			   'count_alleles':count_alleles,
			   'alt_allele':alt_alleles[0] if count_alleles > 1 else 'NA',
			   'depth':sum(counts),
			   'count_ref':sum(counts[0:2]),
			   'count_alt':sum(counts[2:4]),
			   'cons_allele':cons_allele,
			   'ref_freq':'NA'if sum(counts) == 0 else sum(counts[0:2])/float(sum(counts))
			   }

def compute_perc_id(aln):
	""" Compute percent identity for paired-end read """
	length = aln.query_length
	edit = dict(aln.tags)['NM']
	return 100 * (length - edit)/float(length)

## Deprecated
def fetch_pe_reads(aln_file):
	""" Use pysam to yield paired end reads from bam file """
	pe_read = []
	for aln in aln_file.fetch(until_eof = True):
		if not aln.is_paired:
			yield [aln]
		elif aln.mate_is_unmapped and aln.is_read1:
			yield [aln]
		elif aln.mate_is_unmapped and aln.is_read2:
			yield [aln]
		else:
			pe_read.append(aln)
			if len(pe_read) == 2:
				yield pe_read
				pe_read = []


def write_best_hits(args, genome_clusters, best_hits, reference_map, batch_index):
	""" Write reassigned PE reads to disk """
	if args['verbose']: print("    writing mapped reads to disk")
	try: os.makedirs('/'.join([args['out'], 'reassigned']))
	except: pass
	# open filehandles
	aln_files = {}
	scaffold_to_genome = {}
	# loop over genome clusters
	for cluster_id in genome_clusters:
		# get template bam file
		bam_path = '/'.join([args['out'], 'bam', '%s.%s.bam' % (cluster_id, batch_index)])
		if not os.path.isfile(bam_path):
			sys.stderr.write("    bam file not found for %s.%s Skipping\n" % (cluster_id, batch_index))
			continue
		template = pysam.AlignmentFile(bam_path, 'rb')
		# store filehandle
		outpath = '/'.join([args['out'], 'reassigned', '%s.%s.bam' % (cluster_id, batch_index)])
		aln_files[cluster_id] = pysam.AlignmentFile(outpath, 'wb', template=template)
	# write reads to disk
	for cluster_id, pe_read in best_hits.values():
		for aln in pe_read:
			aln_files[cluster_id].write(aln)

def map_reads(args, genome_clusters, batch_index):
	""" find and write bets-hits to disk """
	# Get best hit for each read in batch_index
	best_hits, reference_map = find_best_hits(args, genome_clusters, batch_index)
	# Write best hits to disk
	write_best_hits(args, genome_clusters, best_hits, reference_map, batch_index)

def write_pangene_coverage(args, pangene_to_cov, phyeco_cov, cluster_id):
	""" Write coverage of pangenes for genome cluster to disk """
	outdir = '/'.join([args['out'], 'coverage'])
	try: os.mkdir(outdir)
	except: pass
	outfile = gzip.open('/'.join([outdir, '%s.cov.gz' % cluster_id]), 'w')
	for pangene in sorted(pangene_to_cov.keys()):
		cov = pangene_to_cov[pangene]
		cn = cov/phyeco_cov if phyeco_cov > 0 else 0
		outfile.write('\t'.join([pangene, str(cov), str(cn)])+'\n')

def count_mapped_bp(args):
	""" Count number of bp mapped to each centroid across pangenomes """
	bam_path = '/'.join([args['out'], 'pangenome.bam'])
	aln_file = pysam.AlignmentFile(bam_path, "rb")
	ref_to_length = dict([(i,j) for i,j in zip(aln_file.references, aln_file.lengths)])
	ref_to_cov = dict([(i,0) for i in aln_file.references])
	for aln in aln_file.fetch(until_eof = True):
		query = aln.query_name
		pid = compute_perc_id(aln)
		if pid < args['pid']:
			continue
		else:
			ref_id = aln_file.getrname(aln.reference_id)
			cov = float(aln.query_alignment_length)/ref_to_length[ref_id]
			ref_to_cov[ref_id] += cov
	return ref_to_cov

def compute_phyeco_cov(args, pangene_to_cov, cluster_id):
	""" Compute coverage of phyeco markers for genome cluster """
	markers = ['B000039','B000041','B000062','B000063','B000065','B000071','B000079',
			   'B000080','B000081','B000082','B000086','B000096','B000103','B000114']
	phyeco_covs = []
	inpath = '/'.join([args['db_dir'], cluster_id, 'pangene_to_phyeco.gz'])
	infile = gzip.open(inpath)
	next(infile)
	for line in infile:
		pangene, phyeco_id = line.rstrip().split()
		if phyeco_id in markers:
			phyeco_covs.append(pangene_to_cov[pangene])
	return np.median(phyeco_covs)

def compute_pangenome_coverage(args):
	""" Compute coverage of pangenome for cluster_id and write results to disk """
	ref_to_cov = count_mapped_bp(args)
	phyeco_cov = compute_phyeco_cov(args, pangene_to_cov, cluster_id)
	write_pangene_coverage(args, pangene_to_cov, phyeco_cov, cluster_id)

def max_mem_usage():
	""" Return max mem usage (Gb) of self and child processes """
	max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
	return round((max_mem_self + max_mem_child)/float(1e6), 2)

def convert_to_ascii_quality(scores):
	""" Convert quality scores to Sanger encoded (Phred+33) ascii values """
	ascii = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQR"""
	score_to_ascii = dict((x,y) for x,y in zip(range(0,50),list(ascii)))
	return ''.join([score_to_ascii[x] for x in scores])

def convert_from_ascii_quality(asciis):
	""" Convert quality scores to Sanger encoded (Phred+33) ascii values """
	ascii = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQR"""
	ascii_to_score = dict((y,x) for x,y in zip(range(0,50),list(ascii)))
	return [ascii_to_score[x] for x in asciis]

def write_fastq_record(aln, index, outfile):
	""" Write pysam alignment record to outfile in FASTQ format """
	outfile.write('@%s.%s length=%s\n' % (aln.query_name,str(index),str(aln.query_length)))
	outfile.write('%s\n' % (aln.query_sequence))
	outfile.write('+%s.%s length=%s\n' % (aln.query_name,str(index),str(aln.query_length)))
	outfile.write('%s\n' % convert_to_ascii_quality(aln.query_qualities))

def bam_to_fastq(genome_clusters, args):
	""" Converts bam to fastq for reads assigned to each genome cluster """
	bam_dir = '/'.join([args['out'], 'reassigned'])
	batch_indexes = sorted(set([_.split('.')[1] for _ in os.listdir(bam_dir)]))
	fastq_dir = '/'.join([args['out'], 'fastq'])
	try: os.mkdir(fastq_dir)
	except: pass
	for genome_cluster in genome_clusters:
		outfile = gzip.open(os.path.join(fastq_dir, genome_cluster+'.fastq.gz'), 'w')
		for batch_index in batch_indexes:
			bam_name = '.'.join([genome_cluster, batch_index, 'bam'])
			bam_path = '/'.join([bam_dir, bam_name])
			aln_file = pysam.AlignmentFile(bam_path, "rb")
			for index, aln in enumerate(aln_file.fetch(until_eof = True)):
				write_fastq_record(aln, index, outfile)

def build_pangenome_db(args, genome_clusters):
	""" Build FASTA and BT2 database from pangene cluster centroids """
	# fasta database
	outdir = '/'.join([args['out'], 'db'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	pgdb_fasta = '/'.join([args['out'], 'db', 'pangenomes.fa'])
	outfile = open(pgdb_fasta, 'w')
	indir = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation/pan_genomics2/pan_genomes/90'
	for cluster_id in genome_clusters:
		for line in open('/'.join([indir, cluster_id, 'centroids.fa'])):
			outfile.write(line)
	# bowtie2 database
	pgdb_bt2 = '/'.join([args['out'], 'db', 'pangenomes'])
	command = ' '.join([args['bowtie2-build'], pgdb_fasta, pgdb_bt2])
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()

def build_genome_db(args, genome_clusters):
	""" Build FASTA and BT2 database from genome cluster centroids """
	# fasta database
	outdir = '/'.join([args['out'], 'db'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	genomes_fasta = open('/'.join([args['out'], 'db', 'genomes.fa']), 'w')
	genomes_map = open('/'.join([args['out'], 'db', 'genomes.map']), 'w')
	indir = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation/phylo_cnv/phylo_db/genome_clusters'
	for cluster_id in genome_clusters:
		for line in open('/'.join([indir, cluster_id, 'cluster_centroid/centroid.fna'])):
			genomes_fasta.write(line)
			if line[0] == '>':
				sid = line.rstrip().lstrip('>').split()[0]
				genomes_map.write(sid+'\t'+cluster_id+'\n')
	# bowtie2 database
	inpath = '/'.join([args['out'], 'db', 'genomes.fa'])
	outpath = '/'.join([args['out'], 'db', 'genomes'])
	command = ' '.join([args['bowtie2-build'], inpath, outpath])
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()

def run_pipeline(args):
	""" Run entire pipeline """
	
	check_arguments(args) # Check validity of command line arguments
	add_binaries(args) # Add path to external binaries
	if args['verbose']: print_copyright()

	# 1a. Estimate the abundance of genome-clusters
	if args['species_profile']:
		start = time.time()
		if args['verbose']: print("\nEstimating the abundance of genome-clusters")
		cluster_abundance, cluster_summary = phylo_species.estimate_species_abundance(
			{'inpath':args['m1'], 'nreads':args['reads_ms'],
			 'outpath':'/'.join([args['out'], 'genome_clusters']),
			 'min_quality': 25, 'min_length': 50, 'max_n':0.05,
			 'threads':args['threads']})
		phylo_species.write_abundance('%s/genome_clusters.abundance' % args['out'], cluster_abundance)
		phylo_species.write_summary('%s/genome_clusters.summary' % args['out'], cluster_summary)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 1b. Select genome-clusters for downstream steps
	if args['verbose']: print("\nSelecting genome-clusters for pangenome alignment")
	cluster_abundance = read_phylo_species('/'.join([args['out'], 'genome_clusters.abundance']))
	genome_clusters = select_genome_clusters(cluster_abundance, args)
	if len(genome_clusters) == 0:
		sys.exit("No genome-clusters were detected")
	elif args['verbose']:
		for cluster, abundance in sorted(genome_clusters.items(), key=operator.itemgetter(1), reverse=True):
			print("  cluster_id: %s abundance: %s" % (cluster, round(abundance,2)))

	# 2a. Build pangenome database for selected GCs
	if args['pangenome_build_db']:
		if args['verbose']: print("\nBuilding pangenome database")
		start = time.time()
		build_pangenome_db(args, genome_clusters)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 2a. Use bowtie2 to align reads to pangenome database
	if args['pangenome_align']:
		start = time.time()
		if args['verbose']: print("\nAligning reads to pangenomes")
		pangenome_align(args, args['tax_mask'])
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 2b. Compute pangenome coverage for each genome-cluster
	if args['pangenome_cov']:
		start = time.time()
		if args['verbose']: print("\nComputing coverage of pangenomes")
		compute_pangenome_coverage(args)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 3a. Build genome database for selected GCs
	if args['snps_build_db']:
		if args['verbose']: print("\nBuilding database of representative genomes")
		start = time.time()
		build_genome_db(args, genome_clusters)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 3b. Use bowtie2 to map reads to a representative genome for each genome-cluster
	if args['snps_align']:
		if args['verbose']: print("\nMapping reads to representative genomes")
		start = time.time()
		genome_align(args)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 3c. Use mpileup to identify SNPs
	if args['snps_call']:
		start = time.time()
		if args['verbose']: print("\nRunning mpileup")
		pileup(args)
		format_vcf(args, genome_clusters)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

