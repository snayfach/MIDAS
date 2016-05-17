#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip, numpy as np, random
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import utility, analyze_snps as st

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
Usage: snp_prevalence.py [options]
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
	args = vars(parser.parse_args())
	args.update(st.init_paths(args))
	return args

def compute_prevalence(args, samples):
	snps = set([])
	for sample in samples.values():
		snps |= sample.snps
	prev = dict([(_,0) for _ in snps])
	for sample in samples.values():
		for snp in sample.snps:
			prev[snp] += 1
	write_prevalence(args, prev, samples)

def write_prevalence(args, prevalence, samples):
	outfile = open(args['out'], 'w')
	outfile.write('snp_id\tcount_samples\tfraction_samples\n')
	counts = prevalence.items()
	counts.sort(key=lambda x: x[1], reverse=True)
	for snp_id, count in counts:
		fraction = float(count)/len(samples)
		r = [snp_id, count, fraction]
		outfile.write('\t'.join([str(x) for x in r])+'\n')

if __name__ == '__main__':
	args = parse_arguments()
	utility.print_copyright()
	
	print("Loading samples...")
	samples = st.init_samples(args)
		
	print("Identifying SNPs...")
	st.identify_snps(args, samples)
	
	print("Computing SNP prevalence...")
	prevalence = compute_prevalence(args, samples)








