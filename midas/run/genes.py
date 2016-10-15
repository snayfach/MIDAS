#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, subprocess, gzip
from time import time
from midas import utility
from midas.run import stream_bam

def build_pangenome_db(args, species):
	""" Build FASTA and BT2 database from pangene species centroids """
	import Bio.SeqIO
	# fasta database
	outdir = '/'.join([args['outdir'], 'genes/temp'])
	pangenome_fasta = open('/'.join([outdir, 'pangenomes.fa']), 'w')
	pangenome_map = open('/'.join([outdir, 'pangenome.map']), 'w')
	db_stats = {'total_length':0, 'total_seqs':0, 'species':0}
	for sp in species:
		db_stats['species'] += 1
		infile = utility.iopen(sp.pan_genome)
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
	bampath = '/'.join([args['outdir'], 'genes/temp/pangenome.bam'])
	command += '| %s view -b - > %s' % (args['samtools'], bampath)
	# Run command
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# Check for errors
	print("  finished aligning")
	utility.check_exit_code(process, command)
	print("  checking bamfile integrity")
	utility.check_bamfile(args, bampath)

def count_mapped_bp(args):
	""" Count number of bp mapped to each centroid across pangenomes """
	import pysam, numpy as np
	bam_path = '/'.join([args['outdir'], 'genes/temp/pangenome.bam'])
	aln_file = pysam.AlignmentFile(bam_path, "rb")
	ref_to_length = dict([(i,j) for i,j in zip(aln_file.references, aln_file.lengths)])
	gene_to_cov = dict([(i,0.0) for i in aln_file.references])
	for index, aln in enumerate(aln_file.fetch(until_eof = True)):
		query = aln.query_name
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
			cov = len(aln.query_alignment_sequence)/float(ref_to_length[gene_id])
			gene_to_cov[gene_id] += cov
	return gene_to_cov

def compute_marker_cov(args, species, gene_to_cov, ref_to_species):
	""" Count number of bp mapped to each marker marker gene """
	from numpy import median
	# read in map of gene to marker
	gene_to_marker = read_marker_map(args, species)
	marker_ids = set([marker_id for gene_id, marker_id in gene_to_marker.items()])
	# init marker coverage
	species_to_marker_to_cov = {}
	for species_id in species:
		species_to_marker_to_cov[species_id] = {}
		for marker_id in marker_ids:
			species_to_marker_to_cov[species_id][marker_id] = 0.0
	# compute marker coverages
	for gene_id, marker_id in gene_to_marker.items():
		species_id = ref_to_species[gene_id]
		if marker_id in marker_ids and gene_id in gene_to_cov:
			species_to_marker_to_cov[species_id][marker_id] += gene_to_cov[gene_id]
	# compute median marker cov
	species_to_norm = {}
	for species_id in species_to_marker_to_cov:
		covs = list(species_to_marker_to_cov[species_id].values())
		species_to_norm[species_id] = median(covs)
	return species_to_norm

def compute_pangenome_coverage(args):
	""" Compute coverage of pangenome for species_id and write results to disk """
	
	# map gene_id to species_id
	ref_to_species = {}
	for line in open('/'.join([args['outdir'], 'genes/temp/pangenome.map'])):
		gene_id, species_id = line.rstrip().split()
		ref_to_species[gene_id] = species_id
		
	# open outfiles for each species_id
	outfiles = {}
	species = set(ref_to_species.values())
	for species_id in species:
		outpath = '/'.join([args['outdir'], 'genes/output/%s.genes.gz' % species_id])
		outfiles[species_id] = utility.iopen(outpath, 'w')
		outfiles[species_id].write('\t'.join(['gene_id', 'coverage', 'copy_number'])+'\n')
	
	# parse bam into cov files for each species_id
	gene_to_cov = count_mapped_bp(args)

	# compute normalization factor
	species_to_norm = compute_marker_cov(args, species, gene_to_cov, ref_to_species)

	# write to output files
	for gene_id in sorted(gene_to_cov):
		cov = gene_to_cov[gene_id]
		species_id = ref_to_species[gene_id]
		outfile = outfiles[species_id]
		normcov = cov/species_to_norm[species_id] if species_to_norm[species_id] > 0 else 0.0
		outfile.write('\t'.join([str(x) for x in [gene_id, cov, normcov]])+'\n')

def remove_tmp(args):
	""" Remove specified temporary files """
	import shutil
	shutil.rmtree('/'.join([args['outdir'], 'genes/temp']))

def genes_summary(args):
	""" Get summary of mapping statistics """
	# store stats
	stats = {}
	inpath = '%s/%s' % (args['outdir'], 'genes/temp/pangenome.map')
	for species_id in args['species_id']:
		pangenome_size, covered_genes, total_coverage, marker_coverage = [0,0,0,0]
		for r in utility.parse_file('/'.join([args['outdir'], 'genes/output/%s.genes.gz' % species_id])):
			pangenome_size += 1
			coverage = float(r['coverage'])
			normcov = float(r['copy_number'])
			if coverage > 0:
				covered_genes += 1
				total_coverage += coverage
			if normcov > 0:
				marker_coverage = coverage/normcov
		stats[species_id] = {'pangenome_size':pangenome_size,
							 'covered_genes':covered_genes,
							 'fraction_covered':covered_genes/float(pangenome_size),
							 'mean_coverage':total_coverage/covered_genes if covered_genes > 0 else 0.0,
							 'marker_coverage':marker_coverage}
	# write stats
	fields = ['pangenome_size', 'covered_genes', 'fraction_covered', 'mean_coverage', 'marker_coverage']
	outfile = open('/'.join([args['outdir'], 'genes/summary.txt']), 'w')
	outfile.write('\t'.join(['species_id'] + fields)+'\n')
	for species_id in stats:
		record = [species_id] + [str(stats[species_id][field]) for field in fields]
		outfile.write('\t'.join(record)+'\n')

def read_marker_map(args, species):
	gene_to_marker = {}
	path = '%s/marker_genes/phyeco.map' % args['db']
	for r in utility.parse_file(path):
		if r['species_id'] in species:
			gene_to_marker[r['gene_id']] = r['marker_id']
	return gene_to_marker

class Species:
	""" Base class for species """
	def __init__(self, id):
		self.id = id
		
	def init_ref_db(self, ref_db):
		for ext in ['', '.gz']:
			inpath = '%s/pan_genomes/%s/centroids.ffn%s' % (ref_db, self.id, ext)
			if os.path.isfile(inpath): self.pan_genome = inpath

class Pangenome:
	""" Base class for pan genome """
	def __init__(self):
		pass

def initialize_species(args):
	species = []
	splist = '%s/genes/species.txt' % args['outdir']
	if args['build_db']:
		from midas.run.species import select_species
		with open(splist, 'w') as outfile:
			for id in select_species(args):
				species.append(Species(id))
				outfile.write(id+'\n')
	elif os.path.isfile(splist):
		for line in open(splist):
			species.append(Species(line.rstrip()))
	for sp in species:
		sp.init_ref_db(args['db'])
	return species

def run_pipeline(args):
	""" Run entire pipeline """
	
	# Initialize species
	species = initialize_species(args)
	
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
		compute_pangenome_coverage(args)
		genes_summary(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Optionally remove temporary files
	if args['remove_temp']: remove_tmp(args)
		


