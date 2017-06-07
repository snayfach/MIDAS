#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, subprocess, shutil, csv
import Bio.SeqIO, pysam, numpy as np
from time import time
from midas import utility

class Species:
	""" Base class for species """
	def __init__(self, id):
		self.id = id
		self.paths = {}
		self.aligned_reads = 0
		self.mapped_reads = 0
		self.genome_length = 0
		self.covered_bases = 0
		self.total_depth = 0
		self.fraction_covered = 0
		self.mean_coverage = 0
		
	def fetch_paths(self, ref_db):
		indir =  '%s/rep_genomes/%s' % (ref_db, self.id)
		for ext in ['', '.gz']:
			for type in ['fna', 'features']:
				path = '%s/genome.%s%s' % (indir, type, ext)
				if os.path.isfile(path):
					self.paths[type] = path

class Contig:
	""" Base class for contig """
	def __init__(self, id):
		self.id = id
		
def initialize_species(args):
	species = {}
	splist = '%s/snps/species.txt' % args['outdir']
	if args['build_db']:
		from midas.run.species import select_species
		with open(splist, 'w') as outfile:
			for id in select_species(args):
				species[id] = Species(id)
				outfile.write(id+'\n')
	elif os.path.isfile(splist):
		for line in open(splist):
			id = line.rstrip()
			species[id] = Species(id)
	for sp in species.values():
		sp.fetch_paths(ref_db=args['db'])
	return species

def initialize_contigs(species):
	contigs = {}
	for sp in species.values():
		infile = utility.iopen(sp.paths['fna'])
		for rec in Bio.SeqIO.parse(infile, 'fasta'):
			contig = Contig(rec.id)
			contig.id = rec.id
			contig.seq = str(rec.seq).upper()
			contig.length = len(contig.seq)
			contig.species_id = sp.id
			contigs[contig.id] = contig
		infile.close()
	return contigs
	
def build_genome_db(args, species):
	""" Build FASTA and BT2 database of representative genomes """
	import Bio.SeqIO
	# fasta database
	outfile = open('/'.join([args['outdir'], 'snps/temp/genomes.fa']), 'w')
	db_stats = {'total_length':0, 'total_seqs':0, 'species':0}
	for sp in species.values():
		db_stats['species'] += 1
		infile = utility.iopen(sp.paths['fna'])
		for r in Bio.SeqIO.parse(infile, 'fasta'):
			outfile.write('>%s\n%s\n' % (r.id, str(r.seq).upper()))
			db_stats['total_length'] += len(r.seq)
			db_stats['total_seqs'] += 1
		infile.close()
	outfile.close()
	# print out database stats
	print("  total genomes: %s" % db_stats['species'])
	print("  total contigs: %s" % db_stats['total_seqs'])
	print("  total base-pairs: %s" % db_stats['total_length'])
	# bowtie2 database
	command = '%s ' % args['bowtie2-build']
	command += '--threads %s ' % args['threads']
	command += '%s/snps/temp/genomes.fa ' % args['outdir']
	command += '%s/snps/temp/genomes ' % args['outdir']
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)

def genome_align(args):
	""" Use Bowtie2 to map reads to representative genomes """
	# Bowtie2
	bam_path = os.path.join(args['outdir'], 'snps/temp/genomes.bam')
	command = '%s --no-unal ' % args['bowtie2']
	command += '-x %s ' % '/'.join([args['outdir'], 'snps/temp/genomes']) # index
	if args['max_reads']: command += '-u %s ' % args['max_reads'] # max num of reads
	if args['trim']: command += '--trim3 %s ' % args['trim'] # trim 3'
	command += '--%s ' % args['speed'] # speed/sensitivity
	command += '--threads %s ' % args['threads'] 
	command += '-f ' if args['file_type'] == 'fasta' else '-q ' # input type
	if args['m2']: # -1 and -2 contain paired reads
		command += '-1 %s -2 %s ' % (args['m1'], args['m2']) 
	elif args['interleaved']: # -1 contains paired reads
		command += '--interleaved %s ' % args['m1'] 
	else: # -1 contains unpaired reads
		command += '-U %s ' % args['m1'] 
	# Pipe to samtools
	command += '| %s view -b - ' % args['samtools'] # convert to bam
	command += '--threads %s ' % args['threads']
	command += '| %s sort - ' % args['samtools']
	command += '--threads %s ' % args['threads']
	command += '-o %s ' % bam_path
	# Run command
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# Check for errors
	utility.check_exit_code(process, command)
	print("  finished aligning")
	print("  checking bamfile integrity")
	utility.check_bamfile(args, bam_path)

def index_bam(args):
	start = time()
	print("\nIndexing bamfile")
	args['log'].write("\nIndexing bamfile\n")
	command = '%s index %s/snps/temp/genomes.bam' % (args['samtools'], args['outdir'])
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)
	print("  %s minutes" % round((time() - start)/60, 2) )
	print("  %s Gb maximum memory" % utility.max_mem_usage())

def keep_read(aln):
	global sp
	sp.aligned_reads += 1
	align_len = len(aln.query_alignment_sequence)
	query_len = aln.query_length
	# min pid
	if 100*(align_len-dict(aln.tags)['NM'])/float(align_len) < min_pid:
		return False
	# min read quality
	elif np.mean(aln.query_qualities) < min_readq:
		return False
	# min map quality
	elif aln.mapping_quality < min_mapq:
		return False
	# min aln cov
	elif align_len/float(query_len)  < min_aln_cov:
		return False
	else:
		sp.mapped_reads += 1
		return True

def pysam_pileup(args, species, contigs):
	start = time()
	print("\nCounting alleles")
	args['log'].write("\nCounting alleles\n")
	
	# Set global variables for read filtering
	global sp
	global min_pid
	min_pid = args['mapid']
	global min_readq
	min_readq = args['readq']
	global min_mapq
	min_mapq = args['mapq']
	global min_aln_cov
	min_aln_cov = 0.70
	global min_baseq
	min_baseq = args['baseq']
	
	# open outfiles
	for sp in species.values():
		sp.out = utility.iopen('/'.join([args['outdir'], 'snps/output/%s.snps.gz' % sp.id]), 'w')
		header = ['ref_id', 'ref_pos', 'ref_allele', 'depth', 'count_a', 'count_c', 'count_g', 'count_t']
		sp.out.write('\t'.join(header)+'\n')

	# compute coverage
	bampath = '%s/snps/temp/genomes.bam' % args['outdir']
	with pysam.AlignmentFile(bampath, 'rb') as bamfile:
		for contig in contigs.values():
			
			sp = species[contig.species_id]
			
			counts = bamfile.count_coverage(
				contig.id, 
				start=0, 
				end=contig.length, 
				quality_threshold=min_baseq, 
				read_callback=keep_read)
				
			for i in range(0, contig.length):
				ref_pos = i+1
				ref_allele = contig.seq[i]
				depth = sum([counts[_][i] for _ in range(4)])
				count_a = counts[0][i]
				count_c = counts[1][i]
				count_g = counts[2][i]
				count_t = counts[3][i]
				row = [contig.id, ref_pos, ref_allele, depth, count_a, count_c, count_g, count_t]
				sp.out.write('\t'.join([str(_) for _ in row])+'\n')
				sp.genome_length += 1
				sp.total_depth += depth
				if depth > 0: sp.covered_bases += 1
	
	print("  total aligned reads: %s" % sum([sp.aligned_reads for sp in species.values()]))
	print("  total mapped reads: %s" % sum([sp.mapped_reads for sp in species.values()]))
	
	# coverage summary
	for sp in species.values():
		if sp.genome_length > 0:
			sp.fraction_covered = sp.covered_bases/float(sp.genome_length) 
		if sp.covered_bases > 0:
			sp.mean_coverage = sp.total_depth/float(sp.covered_bases) 

	# close outfiles
	for sp in species.values():
		sp.out.close()

	print("  %s minutes" % round((time() - start)/60, 2) )
	print("  %s Gb maximum memory" % utility.max_mem_usage())

def snps_summary(args, species):
	""" Get summary of mapping statistics """
	
	fields = ['species_id', 'genome_length', 'covered_bases', 'fraction_covered', 'mean_coverage', 'aligned_reads', 'mapped_reads']
	outfile = open(args['outdir'] + '/snps/summary.txt', 'w')
	outfile.write('\t'.join(fields)+'\n')
	
	for sp in species.values():
		outfile.write(sp.id+'\t')
		outfile.write(str(sp.genome_length)+'\t')
		outfile.write(str(sp.covered_bases)+'\t')
		outfile.write(str(sp.fraction_covered)+'\t')
		outfile.write(str(sp.mean_coverage)+'\t')
		outfile.write(str(sp.aligned_reads)+'\t')
		outfile.write(str(sp.mapped_reads)+'\n')
	outfile.close()

def remove_tmp(args):
	""" Remove specified temporary files """
	shutil.rmtree('/'.join([args['outdir'], 'snps/temp']))

def run_pipeline(args):
	""" Run entire pipeline """
		
	# Initialize reference data
	print("\nReading reference data")
	start = time()
	species = initialize_species(args)
	contigs = initialize_contigs(species)
	print("  %s minutes" % round((time() - start)/60, 2) )
	print("  %s Gb maximum memory" % utility.max_mem_usage())
	
	# Build genome database for selected species
	if args['build_db']:
		print("\nBuilding database of representative genomes")
		args['log'].write("\nBuilding database of representative genomes\n")
		start = time()
		build_genome_db(args, species)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Use bowtie2 to map reads to a representative genome for each species
	if args['align']:
		args['file_type'] = utility.auto_detect_file_type(args['m1'])
		print("\nMapping reads to representative genomes")
		args['log'].write("\nMapping reads to representative genomes\n")
		start = time()
		genome_align(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Use mpileup to identify SNPs
	if args['call']:
		index_bam(args)
		pysam_pileup(args, species, contigs)
		snps_summary(args, species)

	# Optionally remove temporary files
	if args['remove_temp']: remove_tmp(args)

