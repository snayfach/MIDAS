#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip, numpy as np, random, analyze_snps as st
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import utility

class Sample:
	""" Base class for samples """
	def __init__(self, id):
		self.id = id
		self.snps = set([])
		self.sites = 0

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: count_snps.py [options]
""")
	parser.add_argument('--out', metavar='PATH', type=str,
		help="""path to output file""")
	parser.add_argument('--inbase', metavar='PATH', type=str, required=True,
		help="""basename for input files""")
	parser.add_argument('--max_sites', metavar='INT', type=int, required=False,
		help="""maximum number of sites to use from input matrix
useful for quick tests (use all)""")
	parser.add_argument('--maf', type=float, metavar='FLOAT', required=False, default=0.05,
		help="""minor allele frequency for SNP definition (0.05)""")
	parser.add_argument('--samples', type=str, metavar='STR', required=False,
		help="""comma-separated list of samples to use (use all)""")
	parser.add_argument('--iter', type=int, metavar='INT', required=False, default=20,
		help="""number of iterations to perform at each sampling depth (20)""")

	args = vars(parser.parse_args())
	args.update(st.init_paths(args))
	return args

def count_snps(args, samples):
	counts = {}
	for i in range(1, len(samples)+1):
		counts[i] = {'union':[], 'xsect':[]}
		for j in range(args['iter']):
			union = set([])
			s = random.sample(samples.values(), i)
			for sample in s:
				union |= sample.snps
			xsect = union.copy()
			for sample in s:
				xsect &= sample.snps
			counts[i]['union'].append(len(union))
			counts[i]['xsect'].append(len(xsect))
	return counts

def write_counts(args, counts):
	outfile = open(args['out'], 'w')
	outfile.write('count_samples\titer\tcount_union\tcount_xsect\n')
	for n in counts:
		union = counts[n]['union']
		xsect = counts[n]['xsect']
		for i, j, k in zip(range(args['iter']), union, xsect):
			r = [n, i, j, k]
			outfile.write('\t'.join([str(x) for x in r])+'\n')


if __name__ == '__main__':
	args = parse_arguments()
	utility.print_copyright()
	
	print("Loading samples...")
	samples = st.init_samples(args)
		
	print("Identifying SNPs...")
	st.identify_snps(args, samples)

	print("Performing SNP rarefaction...")
	counts = count_snps(args, samples)
	write_counts(args, counts)








