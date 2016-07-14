#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, itertools
from operator import itemgetter
from midas import parse

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: 
  1) identify the presence of marker alleles in samples
  2) compute the sharing of marker alleles between pairs of samples
  
Output fields:
  sample1 - identifier for sample 1
  sample2 - identifier for sample 2
  count1 - number of marker alleles found in sample 1
  count2 - number of marker alleles found in sample 2
  count_both - number of marker alleles found in sample 1 and 2
  count_either - number of marker alleles found in sample 1 or 2

Usage: sharing_marker_alleles.py indir outpath markers [options]
""",
	epilog="""""")

	parser.add_argument('indir', type=str,
		help="""path to input directory
contains files: [snps_alt_allele.txt, snps_info.txt, snps_summary.txt
                 snps_depth.txt, snps_ref_freq.txt]""")
	parser.add_argument('outpath', type=str,
		help="""path to output file""")
	parser.add_argument('markers', type=str,
		help="""path to list of marker alleles output by id_marker_alleles.py""")

	parser.add_argument('--max_sites', type=int, metavar='INT',
		help="""maximum number of sites to process (use all)
useful for quick tests""")
	parser.add_argument('--max_samples', type=int, metavar='INT',
		help="""maximum number of samples to process (use all)
useful for quick tests""")

	presabs = parser.add_argument_group("Options for determining allele presence/absence")

	presabs.add_argument('--min_freq', type=float, metavar='FLOAT', default=0.10,
		help="""minimum frequency (proportion of reads) for determining the presence/absence
of a nucleotide at a genomic site (0.10)""")
	presabs.add_argument('--min_reads', type=int, metavar='INT', default=3,
		help="""minimum number of reads for determining the presence/absence
of a nucleotide at a genomic site (3)""")

	args = vars(parser.parse_args())
	if not os.path.isdir(args['indir']): sys.exit("Specified input directory '%s' does not exist" % args['indir'])
	if not os.path.isfile(args['markers']): sys.exit("Specified input file '%s' does not exist" % args['markers'])
	return args

def quantify_markers(args, paths):
	""" determine if marker present in each sample """
	markers = parse.parse_tsv(args['markers']) # generator for marker alleles file
	marker = fetch_marker(markers) # dictionary for 1st marker allele
	sample_ids = [_.split()[0] for _ in open(paths['summary'])][1:] # list of all sample_ids
	if args['max_samples']: sample_ids = sample_ids[0:args['max_samples']] # get subset of samples for testing
	samples = dict([(_, {}) for _ in sample_ids]) # dictionary of marker alleles found in each sample
	for index, site in enumerate(parse.parse_sites(args, paths)):
		if not index % 100000: print("%s sites processed" % index)
		# skip non-discriminative sites
		if (site.ref_id != marker['ref_id']
				or site.ref_pos < marker['ref_pos']):
			continue
		# determine if marker present in each sample
		for sample_id in sample_ids:
			# skip samples without marker
			depth = int(site.depth[sample_id])
			if depth == 0:
				continue
			elif marker['allele'] == site.info['ref_allele']:
				freq = float(site.ref_freq[sample_id])
			elif marker['allele'] == site.alt_allele[sample_id]:
				freq = 1-float(site.ref_freq[sample_id])
			else:
				continue
			# record marker allele for sample_id
			if (freq >= args['min_freq']
					and round(freq * depth) >= args['min_reads']):
				samples[sample_id][site.id] = freq
		# fetch next marker allele
		marker = fetch_marker(markers)
		if marker is None: break # stop when there are no more markers
	return samples

def init_paths(args):
	""" fetch paths to input files """
	paths = {}
	exts = ['alt_allele', 'info', 'summary', 'depth', 'ref_freq']
	for ext in exts:
		inpath = '%s/snps_%s.txt' % (args['indir'], ext)
		if os.path.isfile(inpath):
			paths[ext] = inpath
		else:
			sys.exit("Input file does not exist: %s" % inpath)
	return paths

def fetch_marker(markers):
	try:
		marker = next(markers)
		marker['ref_id'] = marker['site_id'].rsplit('|', 1)[0]
		marker['ref_pos'] = int(marker['site_id'].rsplit('|', 1)[1])
		return marker
	except StopIteration:
		return None

def allele_sharing(x, y):
	a = len(x)
	b = len(y)
	i = len(set(x.keys()) & set(y.keys()))
	u = len(set(x.keys()) | set(y.keys()))
	return a, b, i, u


if __name__ == "__main__":

	args = parse_arguments() # dictionary of args
	
	paths = init_paths(args) # dictionary of paths
		
	# open output file
	outfile = open(args['outpath'], 'w')
	header = ['sample1', 'sample2', 'count1', 'count2', 'count_both', 'count_either']
	outfile.write('\t'.join(header)+'\n')
	
	# determine marker alleles present in each sample
	print("Determining marker alleles present in each sample")
	samples = quantify_markers(args, paths)
	
	# quantify marker allele sharing
	print("Quantifying sharing of marker alleles between samples")
	for index, pair in enumerate(itertools.combinations(samples, r=2)):
		if not index % 500: print("%s sample pairs processed" % index)
		sample1, sample2 = pair
		count1, count2, count_both, count_either = allele_sharing(samples[sample1], samples[sample2])
		record = [sample1, sample2, count1, count2, count_both, count_either]
		outfile.write('\t'.join([str(_) for _ in record])+'\n')

