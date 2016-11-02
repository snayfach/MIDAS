#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, subprocess, gzip, csv, Bio.SeqIO, numpy as np
from collections import defaultdict
from time import time
from midas import utility
from midas.run import stream_bam

class Species:
	""" Base class for species """
	def __init__(self, id):
		self.id = id
		self.paths = {}
		self.genes = []
		self.pangenome_size = 0
		self.reads = 0
		self.bases = 0.0
		self.depth = 0.0
		self.markers = defaultdict(float)
		
	def init_ref_db(self, ref_db):
		""" Set paths to input files """
		self.dir = '%s/pan_genomes/%s' % (ref_db, self.id)
		for ext in ['', '.gz']:
			for file in ['centroids.ffn', 'cluster_info.txt', 'gene_info.txt']:
				inpath = '%s/%s%s' % (self.dir, file, ext)
				if os.path.isfile(inpath):
					self.paths[file] = inpath

def initialize_species(args):
	""" Initialize Species objects """
	species = {}
	splist = '%s/genes/species.txt' % args['outdir']
	if args['build_db']:
		from midas.run.species import select_species
		with open(splist, 'w') as outfile:
			for id in select_species(args):
				species[id] = Species(id)
				outfile.write(id+'\n')
	elif os.path.isfile(splist):
		for line in open(splist):
			species[id] = Species(line.rstrip())
	for sp in species.values():
		sp.init_ref_db(args['db'])
	return species

class Gene:
	""" Base class for gene """
	def __init__(self, id):
		self.id = id
		self.reads = 0
		self.bases = 0.0
		self.depth = 0.0
		self.length = 0
		self.copies = 0.0
		self.marker_id = None

def initialize_genes(args, species):
	""" Initialize Gene objects """
	genes = {}
	# fetch gene_id, species_id, gene length
	for sp in species.values():
		path = sp.paths['centroids.ffn']
		file = utility.iopen(path)
		for seq in Bio.SeqIO.parse(file, 'fasta'):
			genes[seq.id] = Gene(seq.id)
			genes[seq.id].species_id = sp.id
			genes[seq.id].length = len(seq.seq)
			sp.pangenome_size += 1
		file.close()
	# fetch marker_id
	path = '%s/marker_genes/phyeco.map' % args['db']
	file = utility.iopen(path)
	reader = csv.DictReader(file, delimiter='\t')
	for r in reader:
		if r['gene_id'] in genes:
			genes[r['gene_id']].marker_id=r['marker_id']
	file.close()
	return genes

def build_pangenome_db(args, species):
	""" Build FASTA and BT2 database from pangene species centroids """
	import Bio.SeqIO
	# fasta database
	outdir = '/'.join([args['outdir'], 'genes/temp'])
	pangenome_fasta = open('/'.join([outdir, 'pangenomes.fa']), 'w')
	pangenome_map = open('/'.join([outdir, 'pangenomes.map']), 'w')
	db_stats = {'total_length':0, 'total_seqs':0, 'species':0}
	for sp in species.values():
		db_stats['species'] += 1
		infile = utility.iopen(sp.paths['centroids.ffn'])
		for r in Bio.SeqIO.parse(infile, 'fasta'):
			pangenome_fasta.write('>%s\n%s\n' % (r.id, str(r.seq).upper()))
			pangenome_map.write('%s\t%s\n' % (r.id, sp.id))
			db_stats['total_length'] += len(r.seq)
			db_stats['total_seqs'] += 1
		infile.close()
	pangenome_fasta.close()
	pangenome_map.close()
	# print out database stats
	print("  total species: %s" % db_stats['species'])
	print("  total genes: %s" % db_stats['total_seqs'])
	print("  total base-pairs: %s" % db_stats['total_length'])
	# bowtie2 database
	inpath = '/'.join([outdir, 'pangenomes.fa'])
	outpath = '/'.join([outdir, 'pangenomes'])
	command = ' '.join([args['bowtie2-build'], inpath, outpath])
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)

def pangenome_align(args):
	""" Use Bowtie2 to map reads to all specified genome species """
	# Build command
	command = '%s --no-unal ' % args['bowtie2']
	#   index
	command += '-x %s ' % '/'.join([args['outdir'], 'genes/temp/pangenomes'])
	#   specify reads
	if args['max_reads']: command += '-u %s ' % args['max_reads']
	#   trim reads
	if args['trim']: command += '--trim3 %s ' % args['trim']
	#   speed/sensitivity
	command += '--%s-local ' % args['speed']
	#   threads
	command += '--threads %s ' % args['threads']
	#   file type
	if args['file_type'] == 'fasta': command += '-f '
	else: command += '-q '
	#   input file
	if (args['m1'] and args['m2']): command += '-1 %s -2 %s ' % (args['m1'], args['m2'])
	else: command += '-U %s ' % args['m1']
	#   output unsorted bam
	bampath = '/'.join([args['outdir'], 'genes/temp/pangenomes.bam'])
	command += '| %s view -b - > %s' % (args['samtools'], bampath)
	# Run command
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# Check for errors
	utility.check_exit_code(process, command)
	print("  finished aligning")
	print("  checking bamfile integrity")
	utility.check_bamfile(args, bampath)

def pangenome_coverage(args, species, genes):
	""" Compute coverage of pangenome for species_id and write results to disk """
	count_mapped_bp(args, species, genes)
	normalize(args, species, genes)
	write_results(args, species, genes)

def count_mapped_bp(args, species, genes):
	""" Count number of bp mapped to each gene across pangenomes """
	import pysam
	bam_path = '/'.join([args['outdir'], 'genes/temp/pangenomes.bam'])
	aln_file = pysam.AlignmentFile(bam_path, "rb")
	i, j = 0,0
	# loop over alignments, sum values per gene
	for index, aln in enumerate(aln_file.fetch(until_eof = True)):
		i += 1
		if stream_bam.compute_perc_id(aln) < args['mapid']:
			continue
		elif stream_bam.compute_aln_cov(aln) < args['aln_cov']:
			continue
		elif np.mean(aln.query_qualities) < args['readq']:
			continue
		elif aln.mapping_quality < args['mapq']:
			continue
		else:
			gene_id = aln_file.getrname(aln.reference_id)
			aln_len = len(aln.query_alignment_sequence)
			gene_len = genes[gene_id].length
			genes[gene_id].reads += 1
			genes[gene_id].bases += aln_len
			genes[gene_id].depth += aln_len/float(gene_len)
			j += 1
	print("  total aligned reads: %s" % i)
	print("  total mapped reads: %s" % j)
	# loop over genes, sum values per species
	for gene in genes.values():
		species[gene.species_id].reads += gene.reads
		species[gene.species_id].bases += gene.bases
		species[gene.species_id].depth += gene.depth
		species[gene.species_id].genes.append(gene.depth)
	# loop over species, compute summaries
	for sp in species.values():
		sp.covered_genes = sum([1 for _ in sp.genes if _ > 0])
		sp.mean_coverage = np.mean([_ for _ in sp.genes if _ > 0])
		sp.fraction_covered = sp.covered_genes/float(sp.pangenome_size)

def normalize(args, species, genes):
	""" Count number of bp mapped to each marker gene """
	# compute marker depth
	for gene in genes.values():
		if gene.marker_id is not None:
			species[gene.species_id].markers[gene.marker_id] += gene.depth
	# compute median marker depth
	for sp in species.values():
		sp.marker_coverage = np.median(sp.markers.values())
	# normalize genes by median marker depth
	for gene in genes.values():
		sp = species[gene.species_id]
		if sp.marker_coverage > 0:
			gene.copies = gene.depth/sp.marker_coverage

def write_results(args, species, genes):
	""" Write results to disk """
	# open outfiles for each species_id
	header = ['gene_id', 'count_reads', 'coverage', 'copy_number']
	for sp in species.values():
		path = '/'.join([args['outdir'], 'genes/output/%s.genes.gz' % sp.id])
		sp.out = utility.iopen(path, 'w')
		sp.out.write('\t'.join(header)+'\n')
	# write to output files
	for gene_id in sorted(genes):
		gene = genes[gene_id]
		sp = species[gene.species_id]
		values = [gene.id, gene.reads, gene.depth, gene.copies]
		sp.out.write('\t'.join([str(_) for _ in values])+'\n')
	# close output files
	for sp in species.values():
		sp.out.close()
	# summary stats
	path = '/'.join([args['outdir'], 'genes/summary.txt'])
	file = open(path, 'w')
	header = ['species_id', 'pangenome_size', 'covered_genes', 'fraction_covered', 'mean_coverage', 'marker_coverage', 'count_reads']
	file.write('\t'.join(header)+'\n')
	for sp in species.values():
		values = [sp.id, sp.pangenome_size, sp.covered_genes, sp.fraction_covered, sp.mean_coverage, sp.marker_coverage, sp.reads]
		file.write('\t'.join([str(_) for _ in values])+'\n')
	file.close()

def remove_tmp(args):
	""" Remove specified temporary files """
	import shutil
	shutil.rmtree('/'.join([args['outdir'], 'genes/temp']))

def run_pipeline(args):
	""" Run entire pipeline """
	
	# Initialize reference data
	print("\nReading reference data")
	start = time()
	species = initialize_species(args)
	genes = initialize_genes(args, species)
	print("  %s minutes" % round((time() - start)/60, 2) )
	print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Build pangenome database for selected species
	if args['build_db']:
		print("\nBuilding pangenome database")
		args['log'].write("\nBuilding pangenome database\n")
		start = time()
		build_pangenome_db(args, species)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Use bowtie2 to align reads to pangenome database
	if args['align']:
		start = time()
		print("\nAligning reads to pangenomes")
		args['log'].write("\nAligning reads to pangenomes\n")
		pangenome_align(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Compute pangenome coverage for each species
	if args['cov']:
		start = time()
		print("\nComputing coverage of pangenomes")
		args['log'].write("\nComputing coverage of pangenomes\n")
		pangenome_coverage(args, species, genes)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Optionally remove temporary files
	if args['remove_temp']: remove_tmp(args)
		


