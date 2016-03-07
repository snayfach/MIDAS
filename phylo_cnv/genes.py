#!/usr/bin/python

# PhyloCNV - estimating the abundance, gene-content, and phylogeny of microbes from metagenomes
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.2'

# Libraries
# ---------
import sys
import os
import subprocess
import gzip
from time import time
import platform
import utility
import filter_bam

# Functions
# ---------

def pangenome_align(args):
	""" Use Bowtie2 to map reads to all specified genome clusters """
	# Build command
	command = '%s --no-unal ' % args['bowtie2']
	#   index
	command += '-x %s ' % '/'.join([args['outdir'], 'db', 'pangenomes'])
	#   specify reads
	if args['reads']: command += '-u %s ' % args['reads']
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
	else: command += '-U %s' % args['m1']
	#   output unsorted bam
	bampath = '/'.join([args['outdir'], 'pangenome.bam'])
	command += '| %s view -b - > %s' % (args['samtools'], bampath)
	# Run command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# Check for errors
	utility.check_exit_code(process, command)
	utility.check_bamfile(args, bampath)

def read_ref_to_cluster(args, type):
	""" Read in map of gene id to genome-cluster id """
	ref_to_cluster = {}
	for line in open('/'.join([args['outdir'], 'db/%s.map' % type])):
		ref_id, cluster_id = line.rstrip().split()
		ref_to_cluster[ref_id] = cluster_id
	return ref_to_cluster

def genes_summary(args):
	""" Get summary of mapping statistics """
	# store stats
	stats = {}
	for cluster_id in set(read_ref_to_cluster(args, 'pangenome').values()):
		pangenome_size, covered_genes, total_coverage, phyeco_coverage = [0,0,0,0]
		for r in utility.parse_file('/'.join([args['outdir'], 'coverage/%s.cov.gz' % cluster_id])):
			pangenome_size += 1
			coverage = float(r['raw_coverage'])
			normcov = float(r['normalized_coverage'])
			if coverage > 0:
				covered_genes += 1
				total_coverage += coverage
			if normcov > 0:
				phyeco_coverage = coverage/normcov
		stats[cluster_id] = {'pangenome_size':pangenome_size,
							 'covered_genes':covered_genes,
							 'mean_coverage':total_coverage/covered_genes if covered_genes > 0 else 0.0,
							 'phyeco_coverage':phyeco_coverage}
	# write stats
	fields = ['pangenome_size', 'covered_genes', 'mean_coverage', 'phyeco_coverage']
	outfile = open('/'.join([args['outdir'], 'genes_summary_stats.txt']), 'w')
	outfile.write('\t'.join(['cluster_id'] + fields)+'\n')
	for cluster_id in stats:
		record = [cluster_id] + [str(stats[cluster_id][field]) for field in fields]
		outfile.write('\t'.join(record)+'\n')

def count_mapped_bp(args):
	""" Count number of bp mapped to each centroid across pangenomes """
	import pysam, numpy as np
	bam_path = '/'.join([args['outdir'], 'pangenome.bam'])
	aln_file = pysam.AlignmentFile(bam_path, "rb")
	ref_to_length = dict([(i,j) for i,j in zip(aln_file.references, aln_file.lengths)])
	ref_to_cov = dict([(i,0.0) for i in aln_file.references])
	for index, aln in enumerate(aln_file.fetch(until_eof = True)):
		query = aln.query_name
		if filter_bam.compute_perc_id(aln) < args['mapid']:
			continue
		elif filter_bam.compute_aln_cov(aln) < args['aln_cov']:
			continue
		elif np.mean(aln.query_qualities) < args['readq']:
			continue
		elif aln.mapping_quality < args['mapq']: # TEST THIS
			continue
		else:
			ref_id = aln_file.getrname(aln.reference_id)
			cov = len(aln.query_alignment_sequence)/float(ref_to_length[ref_id])
			ref_to_cov[ref_id] += cov
	return ref_to_cov

def compute_phyeco_cov(args, genome_clusters, ref_to_cov, ref_to_cluster):
	""" Count number of bp mapped to each PhyEco marker gene """
	from numpy import median
	# read in set of phyeco markers for normalization
	phyeco_ids = set([])
	for line in open(args['pid_cutoffs']):
		phyeco_id, pid = line.rstrip().split()
		phyeco_ids.add(phyeco_id)
	# read in map of gene to phyeco marker
	ref_to_phyeco = {}
	for cluster_id in genome_clusters:
		inpath = '/'.join([args['db'], 'genome_clusters', cluster_id, 'universal_genes.txt.gz'])
		infile = gzip.open(inpath)
		next(infile)
		for line in infile:
			gene_id, phyeco_id = line.rstrip().split()
			ref_to_phyeco[gene_id] = phyeco_id
	# init phyeco coverage
	cluster_to_phyeco_to_cov = {}
	for cluster_id in genome_clusters:
		cluster_to_phyeco_to_cov[cluster_id] = {}
		for phyeco_id in phyeco_ids:
			cluster_to_phyeco_to_cov[cluster_id][phyeco_id] = 0.0
	# compute phyeco coverages
	for ref_id, phyeco_id in ref_to_phyeco.items():
		cluster_id = ref_to_cluster[ref_id]
		if phyeco_id in phyeco_ids and ref_id in ref_to_cov:
			cluster_to_phyeco_to_cov[cluster_id][phyeco_id] += ref_to_cov[ref_id]
	# compute median phyeco cov
	cluster_to_norm = {}
	for cluster_id in cluster_to_phyeco_to_cov:
		covs = cluster_to_phyeco_to_cov[cluster_id].values()
		cluster_to_norm[cluster_id] = median(covs)
	return cluster_to_norm

def compute_pangenome_coverage(args):
	""" Compute coverage of pangenome for cluster_id and write results to disk """
	# map ref_id to cluster_id
	ref_to_cluster = {}
	for line in open('/'.join([args['outdir'], 'db/pangenome.map'])):
		ref_id, cluster_id = line.rstrip().split()
		ref_to_cluster[ref_id] = cluster_id
	# open outfiles for each cluster_id
	outdir = '/'.join([args['outdir'], 'coverage'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	outfiles = {}
	genome_clusters = set(ref_to_cluster.values())
	for cluster_id in genome_clusters:
		outfiles[cluster_id] = gzip.open('/'.join([outdir, '%s.cov.gz' % cluster_id]), 'w')
		outfiles[cluster_id].write('\t'.join(['gene_id', 'raw_coverage', 'normalized_coverage'])+'\n')
	# parse bam into cov files for each cluster_id
	ref_to_cov = count_mapped_bp(args)
	# compute normalization factor
	cluster_to_norm = compute_phyeco_cov(args, genome_clusters, ref_to_cov, ref_to_cluster)
	# write to output files
	for ref_id in sorted(ref_to_cov):
		cov = ref_to_cov[ref_id]
		cluster_id = ref_to_cluster[ref_id]
		outfile = outfiles[cluster_id]
		normcov = cov/cluster_to_norm[cluster_id] if cluster_to_norm[cluster_id] > 0 else 0
		outfile.write('\t'.join([str(x) for x in [ref_id, cov, normcov]])+'\n')

def build_pangenome_db(args, genome_clusters):
	""" Build FASTA and BT2 database from pangene cluster centroids """
	import Bio.SeqIO
	# fasta database
	outdir = '/'.join([args['outdir'], 'db'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	pangenome_fasta = open('/'.join([args['outdir'], 'db/pangenomes.fa']), 'w')
	pangenome_map = open('/'.join([args['outdir'], 'db/pangenome.map']), 'w')
	db_stats = {'total_length':0, 'total_seqs':0, 'genome_clusters':0}
	for cluster_id in genome_clusters:
		db_stats['genome_clusters'] += 1
		inpath = '/'.join([args['db'], 'genome_clusters', cluster_id, 'pangenome.fa.gz'])
		infile = gzip.open(inpath)
		for r in Bio.SeqIO.parse(infile, 'fasta'):
			genome_id = '.'.join(r.id.split('.')[0:2])
			if not args['tax_mask'] or genome_id not in args['tax_mask']:
				pangenome_fasta.write('>%s\n%s\n' % (r.id, str(r.seq)))
				pangenome_map.write('%s\t%s\n' % (r.id, cluster_id))
				db_stats['total_length'] += len(r.seq)
				db_stats['total_seqs'] += 1
	# print out database stats
	print("  total genome-clusters: %s" % db_stats['genome_clusters'])
	print("  total genes: %s" % db_stats['total_seqs'])
	print("  total base-pairs: %s" % db_stats['total_length'])
	# bowtie2 database
	inpath = '/'.join([args['outdir'], 'db/pangenomes.fa'])
	outpath = '/'.join([args['outdir'], 'db/pangenomes'])
	command = ' '.join([args['bowtie2-build'], inpath, outpath])
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)

def remove_tmp(args):
	""" Remove specified temporary files """
	import shutil
	shutil.rmtree('/'.join([args['outdir'], 'db']))
	os.remove('%s/pangenome.bam' % args['outdir'])
		
def run_pipeline(args):
	""" Run entire pipeline """
	
	utility.add_executables(args) # Add paths to external files and binaries
	utility.add_data_files(args)
	
	# Build pangenome database for selected GCs
	if args['build_db']:
		import species
		print("\nBuilding pangenome database")
		start = time()
		genome_clusters = species.select_genome_clusters(args)
		build_pangenome_db(args, genome_clusters)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory") % utility.max_mem_usage()

	# Use bowtie2 to align reads to pangenome database
	if args['align']:
		start = time()
		print("\nAligning reads to pangenomes")
		args['file_type'] = utility.auto_detect_file_type(args['m1'])
		pangenome_align(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory") % utility.max_mem_usage()

	# Compute pangenome coverage for each genome-cluster
	if args['cov']:
		start = time()
		print("\nComputing coverage of pangenomes")
		compute_pangenome_coverage(args)
		genes_summary(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory") % utility.max_mem_usage()

	# Optionally remove temporary files
	if args['remove_temp']: remove_tmp(args)
		


