#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

## TO DO
# exclude genes (?)

__version__ = '0.0.2'

import argparse, sys, os, gzip

def parse_arguments():
	""" Parse command line arguments """
	
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
		
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', '--indir', type=str, dest='in', help='input directory', required=True)
	io.add_argument('-o', '--outdir', type=str, dest='out', help='output directory', required=True)
	io.add_argument('-g', '--genome_cluster', type=str,  help='genome cluster id', required=True)
	
	sample = parser.add_argument_group('Sample filters')
	sample.add_argument('--sample_list', dest='sample_list', type=str,
		default=None, help='file of sample ids to include; each line should contain one id')
	sample.add_argument('--marker_coverage', type=float,
		default=1.0, help='min read depth per sample across 15 phylogenetic marker genes (1.0)')
	sample.add_argument('--gene_coverage', type=float,
		default=0.0, help='min read depth per sample across all genes with non-zero coverage (0.0)')

	gene = parser.add_argument_group('Presence/Absence')
	gene.add_argument('--min_copy', type=float,
		default=0.35, help='genes >= MIN_COPY: classified as present in sample genes < MIN_COPY: classified as absent in sample (0.35)')
	
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
	for sample_id in os.listdir(args['in']):
		# read in summary stats for sample across genome-clusters
		indir = '/'.join([args['in'], sample_id])
		inpath = '/'.join([indir, 'genes_summary_stats.txt'])
		if not os.path.isfile(inpath):
			continue
		else:
			genes_summary = read_file(inpath)
		# check whether sample passes QC
		if not os.path.isfile('%s/coverage/%s.cov.gz' % (indir, args['genome_cluster'])):
			continue
		elif float(genes_summary[args['genome_cluster']]['phyeco_coverage']) < args['marker_coverage']:
			continue
		elif float(genes_summary[args['genome_cluster']]['mean_coverage']) < args['gene_coverage']:
			continue
		# sample passes qc
		else:
			samples.append(sample_id)
	if len(samples) == 0:
		sys.exit("Error: no samples met selection criteria!")
	return samples

def parse_genes(inpath):
	""" Yields formatted records from coverage output """
	infile = gzip.open(inpath)
	fields = next(infile).rstrip().split()
	for line in infile:
		yield dict([(i,j) for i,j in zip(fields, line.rstrip().split())])

def open_infiles(args, samples):
	""" Open coverage files for genome-cluster across all samples """
	infiles = {}
	for sample_id in samples:
		inpath = '%s/%s/coverage/%s.cov.gz' % (args['in'], sample_id, args['genome_cluster'])
		infiles[sample_id] = parse_genes(inpath)
	return infiles

def open_outfiles(args, samples):
	""" Open output files and write headers """
	# open outfiles
	outfiles = {}
	outfiles['coverage_matrix'] = open('%s/%s.coverage' % (args['out'], args['genome_cluster']), 'w')
	outfiles['copynum_matrix'] = open('%s/%s.copynum' % (args['out'], args['genome_cluster']), 'w')
	outfiles['presabs_matrix'] = open('%s/%s.presabs' % (args['out'], args['genome_cluster']), 'w')
	# write headers
	header = ['gene_id'] + samples
	outfiles['coverage_matrix'].write('\t'.join(header)+'\n')
	outfiles['copynum_matrix'].write('\t'.join(header)+'\n')
	outfiles['presabs_matrix'].write('\t'.join(header)+'\n')
	return outfiles

def fetch_gene(infiles, samples):
	""" Fetch SNP data across samples """
	gene = []
	for sample_id in samples:
		gene.append(next(infiles[sample_id]))
	return gene

def presence_absence(gene, min_copynum):
	""" Compute presence/absences from normlized coverage values """
	presabs = ['1' if float(g['normalized_coverage']) >= min_copynum else '0' for g in gene]
	return presabs

def write_records(gene, outfiles, args):
	""" Write snp to outfiles """
	id = [gene[0]['ref_id']]
	outfiles['coverage_matrix'].write('\t'.join(id+[g['raw_coverage'] for g in gene])+'\n')
	outfiles['copynum_matrix'].write('\t'.join(id+[g['normalized_coverage'] for g in gene])+'\n')
	outfiles['presabs_matrix'].write('\t'.join(id+presence_absence(gene, args['min_copy']))+'\n')

def merge_genes(args, infiles, outfiles, samples):
	""" Merge genes across selected samples """
	inpath = '%s/%s/coverage/%s.cov.gz' % (args['in'], infiles.keys()[0], args['genome_cluster'])
	dummyfile = gzip.open(inpath)
	next(dummyfile)
	for line in dummyfile: # outer loop: stop when eof reached
		gene = fetch_gene(infiles, samples) # inner loop: get data for snp across all samples
		write_records(gene, outfiles, args) # write data to output files
	# close open file handles
	for file in infiles.values(): file.close()
	for file in outfiles.values(): file.close()

if __name__ == '__main__':

	args = parse_arguments()
	if not os.path.isdir(args['out']): os.mkdir(args['out'])

	# id samples with sufficient depth
	print("Identifying samples")
	samples = identify_samples(args)
	
	print("Opening input files")
	infiles = open_infiles(args, samples)
	
	print("Opening output files")
	outfiles = open_outfiles(args, samples)

	print("Merging genes and writing output files")
	merge_genes(args, infiles, outfiles, samples)

			
