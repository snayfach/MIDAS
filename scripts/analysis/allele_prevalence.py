#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os
from midas import parse

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="Usage: test.py indir [options]")
	parser.add_argument('indir', type=str,
		help="""path to input directory
contains files: [snps_alt_allele.txt, snps_info.txt, snps_summary.txt
                 snps_depth.txt, snps_ref_freq.txt]""")
	parser.add_argument('outpath', type=str,
		help="""path to output file""")
	parser.add_argument('samples', type=str,
		help="""comma-separated list of sample ids present in INDIR
ex: sample1,sample2,sample3""")
	parser.add_argument('--max_sites', type=int, metavar='INT',
		help="""maximum number of genomic sites to process (use all)
useful for quick tests""")
	parser.add_argument('--min_freq', type=int, metavar='FLOAT', default=0.10,
		help="""minimum frequency for determining the presence/absence 
of a nucleotide at a genomic site (0.10)""")
	args = vars(parser.parse_args())
	args['samples'] = args['samples'].rstrip(',').split(',')
	return args

def allele_props(site, samples, min_freq):
	""" compute the prevelance of the 4 nucleotides for site across samples """
	total = 0
	counts = {'A':0, 'T':0, 'C':0, 'G':0}
	for sample_id in samples: # fetch site info for each sample
		if int(site.depth[sample_id]) == 0: # skip samples with no reads
			continue
		ref_freq = float(site.ref_freq[sample_id]) # count reference allele
		if ref_freq >= min_freq:
			ref_allele = site.info['ref_allele']
			counts[ref_allele] += 1
		alt_freq = 1-ref_freq # count alternate allele
		if alt_freq >= min_freq:
			alt_allele = site.alt_allele[sample_id]
			counts[alt_allele] += 1
		total += 1 # keep track of number of samples with data
	# normalize proportions
	props = {}
	for allele, count in counts.items():
		props[allele] = float(count)/total if total > 0 else 0.0
	return props

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

def check_sample_ids(args, paths):
	""" make sure that specified sample ids are present in input """
	samples = [_.split()[0] for _ in open(paths['summary'])][1:]
	for sample_id in args['samples']:
		if sample_id not in samples:
			sys.exit("Specified sample '%s' not present in indir" % sample_id)

if __name__ == "__main__":

	args = parse_arguments()
	paths = init_paths(args)
	check_sample_ids(args, paths)

	bases = ['A', 'T', 'C', 'G']
	outfile = open(args['outpath'], 'w')
	outfile.write('site_id\t' + '\t'.join(bases) + '\n')

	for site in parse.parse_sites(args, paths):
		props = allele_props(site, args['samples'], args['min_freq'])
		outfile.write(site.id+'\t')
		outfile.write('\t'.join([str(props[base]) for base in bases]) + '\n')





