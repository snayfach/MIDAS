#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

## TO DO
# exclude genes (?)

__version__ = '0.0.2'

import argparse, sys, os, gzip
from collections import defaultdict

def parse_arguments():
	""" Parse command line arguments """
	
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
		
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', '--indir', type=str, dest='in', help='input directory', required=True)
	io.add_argument('-o', '--outdir', type=str, dest='out', help='output directory', required=True)
	io.add_argument('-g', '--genome_cluster', type=str,  help='genome cluster id', required=True)
	io.add_argument('-D', '--db', type=str,  help='directory to genome-clusters database', required=True)
	
	sample = parser.add_argument_group('Sample filters')
	sample.add_argument('--sample_list', dest='sample_list', type=str,
		default=None, help='file of sample ids to include; each line should contain one id')
	sample.add_argument('--marker_coverage', type=float,
		default=1.0, help='min read depth per sample across 15 phylogenetic marker genes (1.0)')
	sample.add_argument('--gene_coverage', type=float,
		default=0.0, help='min read depth per sample across all genes with non-zero coverage (0.0)')
	sample.add_argument('--max_samples', type=int, help='maximum number of samples to process; useful for testing (use all)')

	gene = parser.add_argument_group('Presence/Absence')
	gene.add_argument('--min_copy', type=float,
		default=0.35, help='genes >= MIN_COPY: classified as present in sample genes < MIN_COPY: classified as absent in sample (0.35)')
	gene.add_argument('--cluster_pid', type=str, dest='cluster_pid',
		default='95', choices=['75', '80', '85', '90', '95', '99'],
		help='Gene family percent identity. Small values => fewer, larger gene families. Large values => more, smaller gene families')
	
	args = vars(parser.parse_args())
	
	return args

def read_file(inpath):
	""" Read in summary snp statistics for genome-clusters """
	d = {}
	infile = open(inpath)
	fields = next(infile).rstrip().split()
	for line in open(inpath):
		values = line.rstrip().split()
		rec = dict([(i,j) for i,j in zip(fields, values)])
		d[rec['cluster_id']] = rec
	return d

def identify_samples(args):
	""" Identify samples with sufficient depth for target genome-cluster """
	samples = []
	sample_list = [x.rstrip() for x in open(args['sample_list']).readlines()] if args['sample_list'] else None
	for sample_id in os.listdir(args['in']):
		# read in summary stats for sample across genome-clusters
		indir = '/'.join([args['in'], sample_id])
		inpath = '/'.join([indir, 'genes_summary_stats.txt'])
		if not os.path.isfile(inpath):
			continue
		else:
			genes_summary = read_file(inpath)
		# check whether sample passes QC
		if args['sample_list'] and sample_id not in sample_list:
			continue
		elif args['genome_cluster'] not in genes_summary:
			continue
		elif float(genes_summary[args['genome_cluster']]['phyeco_coverage']) < args['marker_coverage']:
			continue
		elif float(genes_summary[args['genome_cluster']]['mean_coverage']) < args['gene_coverage']:
			continue
		# sample passes qc
		else:
			samples.append(sample_id)
			if args['max_samples'] and len(samples) >= args['max_samples']: # only keep max_samples is specified
				break
	if len(samples) == 0:
		sys.exit("Error: no samples met selection criteria!")
	return samples

def parse_genes(inpath):
	""" Yields formatted records from coverage output """
	infile = gzip.open(inpath)
	fields = next(infile).rstrip().split()
	for line in infile:
		yield dict([(i,j) for i,j in zip(fields, line.rstrip().split())])



def read_fam_map(args):
	""" Read gene family map """
	fam_map = {}
	inpath = '%s/%s/gene_family_map.txt.gz' % (args['db'], args['genome_cluster'])
	infile = gzip.open(inpath)
	fields = next(infile).rstrip().split()
	for line in infile:
		values = line.rstrip().split()
		map = dict([(f,v) for f,v in zip(fields, values)])
		fam_map[map['99']] = map[args['cluster_pid']]
	return fam_map

def compute_sample_copy_num(samples, args, fam_map):
	""" Compute gene copy numbers for samples """
	fam_copy_num = {}
	for sample_id in samples:
		d = defaultdict(float)
		inpath = '%s/%s/coverage/%s.cov.gz' % (args['in'], sample_id, args['genome_cluster'])
		for r in parse_genes(inpath):
			fam_id = fam_map[r['ref_id']]
			d[fam_id] += float(r['normalized_coverage'])
		fam_copy_num[sample_id] = d
	return fam_copy_num

def read_genome_ids(args):
	""" Read in genome ids for genome cluster """
	genome_ids = set([])
	inpath = '%s/%s/genomes.txt.gz' % (args['db'], args['genome_cluster'])
	infile = gzip.open(inpath); next(infile)
	for line in infile:
		genome_id = line.split('\t')[1]
		genome_ids.add(genome_id)
	return list(genome_ids)

def compute_ref_copy_num(genome_ids, args, fam_map):
	""" Compute gene copy numbers for reference genomes """
	ref_copy_num = dict([(_,defaultdict(float)) for _ in genome_ids])
	inpath = '%s/%s/ref_copy_num.txt.gz' % (args['db'], args['genome_cluster'])
	infile = gzip.open(inpath); next(infile)
	for line in infile:
		genome_id, gene_id, copy_number = line.rstrip().split()
		fam_id = fam_map[gene_id]
		ref_copy_num[genome_id][fam_id] += float(copy_number)
	return ref_copy_num

def presence_absence(copy_nums, min_copynum):
	""" Compute presence/absences from normlized coverage values """
	presabs = [1 if copy_num >= min_copynum else 0 for copy_num in copy_nums]
	return presabs

def write_results(samples, genome_ids, sample_copy_num, ref_copy_num, args):
	# open outfiles
	outfiles = {}
	outfiles['copynum'] = open('%s/%s.copynum' % (args['out'], args['genome_cluster']), 'w')
	outfiles['presabs'] = open('%s/%s.presabs' % (args['out'], args['genome_cluster']), 'w')
	# write headers
	header = ['gene_id'] + samples + genome_ids
	outfiles['copynum'].write('\t'.join(header)+'\n')
	outfiles['presabs'].write('\t'.join(header)+'\n')
	# write values
	for gene_id in genes:
		copy_nums = []
		outfiles['copynum'].write('%s\t' % gene_id) # write gene id
		outfiles['presabs'].write('%s\t' % gene_id)
		for sample_id in samples:
			try: copy_nums.append(sample_copy_num[sample_id][gene_id])
			except: copy_nums.append(0.0)
		for genome_id in genome_ids:
			try: copy_nums.append(ref_copy_num[genome_id][gene_id])
			except: copy_nums.append(0.0)
		presabs = presence_absence(copy_nums, args['min_copy'])
		outfiles['copynum'].write('\t'.join([str(_) for _ in copy_nums])+'\n') # write gene values
		outfiles['presabs'].write('\t'.join([str(_) for _ in presabs])+'\n')

if __name__ == '__main__':

	args = parse_arguments()
	if not os.path.isdir(args['out']): os.mkdir(args['out'])

	print("Mapping gene ids")
	fam_map = read_fam_map(args) # map 99% gene ids to lower level
	genes = set(fam_map.values())
	
	print("Computing gene copy numbers for samples")
	samples = identify_samples(args) # id samples with sufficient depth
	sample_copy_num = compute_sample_copy_num(samples, args, fam_map) # gene copy numbers across samples
	
	print("Computing gene copy numbers for reference genomes")
	genome_ids = read_genome_ids(args)
	ref_copy_num = compute_ref_copy_num(genome_ids, args, fam_map) # gene copy numbers across reference genomes

	print("Writing results")
	write_results(samples, genome_ids, sample_copy_num, ref_copy_num, args)



