#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, subprocess, shutil, csv
import Bio.SeqIO, pysam, numpy as np
from time import time
from midas import utility
from smelter.iggdb import IGGdb
from smelter.utilities import tsprint
from multiprocessing import Queue
import multiprocessing
import json

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

	def fetch_paths(self, iggdb):
		self.paths['fna'] = iggdb.get_species(species_id=self.id)['representative_genome_path']

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
	iggdb = args['iggdb']
	for sp in species.values():
		sp.fetch_paths(iggdb)
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
	command += '--%s' % args['speed'] # alignment speed
	command += '-local ' if args['mode'] == 'local' else ' ' # global/local alignment
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

def keep_read_work(aln, my_args, aln_stats):
	aln_stats['aligned_reads'] += 1
	# align and query length
	align_len = len(aln.query_alignment_sequence)
	query_len = aln.query_length
	# min pid filter
	if 100*(align_len-dict(aln.tags)['NM'])/float(align_len) < my_args['mapid']:
		return False
	# min read quality filter
	elif np.mean(aln.query_qualities) < my_args['readq']:
		return False
	# min map quality filter
	elif aln.mapping_quality < my_args['mapq']:
		return False
	# min aln cov filter
	elif align_len/float(query_len)  < my_args['aln_cov']:
		return False
	# read passes all filters
	else:
		aln_stats['mapped_reads'] += 1
		return True


def species_pileup(species_id, contigs):

	global global_args
	args = global_args

	# summary stats
	aln_stats = {'genome_length':0,
				 'total_depth':0,
				 'covered_bases':0,
				 'aligned_reads':0,
				 'mapped_reads':0}

	def keep_read(x):
		return keep_read_work(x, global_args, aln_stats)

	# open outfiles
	out_path = '%s/snps/output/%s.snps.gz' % (args['outdir'], species_id)
	out_file = utility.iopen(out_path, 'w')
	header = ['ref_id', 'ref_pos', 'ref_allele', 'depth', 'count_a', 'count_c', 'count_g', 'count_t']
	out_file.write('\t'.join(header)+'\n')

	# compute coverage
	bampath = '%s/snps/temp/genomes.bam' % args['outdir']
	with pysam.AlignmentFile(bampath, 'rb') as bamfile:
		for contig_id in sorted(list(contigs.keys())):

			contig = contigs[contig_id]

			if contig['species_id'] != species_id:
				continue

			counts = bamfile.count_coverage(
				contig_id,
				start=0,
				end=contig['length'],
				quality_threshold=args['baseq'],
				read_callback=keep_read)

			for i in range(0, contig['length']):
				ref_pos = i+1
				ref_allele = contig['seq'][i]
				depth = sum([counts[_][i] for _ in range(4)])
				count_a = counts[0][i]
				count_c = counts[1][i]
				count_g = counts[2][i]
				count_t = counts[3][i]
				row = [contig_id, ref_pos, ref_allele, depth, count_a, count_c, count_g, count_t]
				out_file.write('\t'.join([str(_) for _ in row])+'\n')
				aln_stats['genome_length'] += 1
				aln_stats['total_depth'] += depth
				if depth > 0: aln_stats['covered_bases'] += 1

	out_file.close()
	tsprint(json.dumps({species_id: aln_stats}, indent=4))
	return (species_id, aln_stats)


def pysam_pileup(args, species, contigs):
	start = time()
	print("\nCounting alleles")
	args['log'].write("\nCounting alleles\n")

	# We cannot pass args to a subprocess unfortunately because args['log'] is an object;
	# so we can make it a global, although that is certainly living dangerously.
	# TODO: Just clean this up.
	global global_args
	global_args = args

	# run pileups per species in parallel
	argument_list = []
	# We might not need this for contigs.  It was an attempt to eliminate the nonserializable subprocess argument.  Which is args.
	contigs = { str(c.id): {'species_id': str(c.species_id), 'length': int(c.length), 'seq': [chr for chr in c.seq]} for c in contigs.values() }
	for species_id in species:
		argument_list.append([species_id, contigs])

	mp = multiprocessing.Pool(int(args['threads']))
	# update alignment stats for species objects
	for species_id, stats in mp.starmap(species_pileup, argument_list):
		sp = species[species_id]
		sp.genome_length = stats['genome_length']
		sp.covered_bases = stats['covered_bases']
		sp.total_depth = stats['total_depth']
		sp.aligned_reads = stats['aligned_reads']
		sp.mapped_reads = stats['mapped_reads']
		if sp.genome_length > 0:
			sp.fraction_covered = sp.covered_bases/float(sp.genome_length)
		if sp.covered_bases > 0:
			sp.mean_coverage = sp.total_depth/float(sp.covered_bases)

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
	if 'db' in args:
		args['iggdb'] = IGGdb(f"{args['db']}/metadata/species_info.tsv")
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
