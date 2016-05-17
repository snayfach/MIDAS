#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip, numpy as np, random
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import utility, analyze_snps as st

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: snp_sharing.py [options]
""")
	parser.add_argument('--out', metavar='PATH', type=str,
		help="""path to output file""")
	parser.add_argument('--inbase', metavar='PATH', type=str, required=True,
		help="""basename for input files""")
	parser.add_argument('--max_sites', metavar='INT', type=int, required=False,
		help="""maximum number of sites to use from input matrix
useful for quick tests (use all)""")
#	parser.add_argument('--min_ratio', type=float, metavar='FLOAT', required=False, default=float('-inf'),
#		help="""minimum ratio of per-site depth versus genome-average (-inf)""")
#	parser.add_argument('--max_ratio', type=float, metavar='FLOAT', required=False, default=float('inf'),
#		help="""maximum ratio of per-site depth versus genome-average (inf)""")
	parser.add_argument('--maf', type=float, metavar='FLOAT', required=False, default=0.05,
		help="""minor allele frequency for SNP definition (0.05)""")
	parser.add_argument('--samples', type=str, metavar='STR', required=False,
		help="""comma-separated list of samples to use (use all)""")
	args = vars(parser.parse_args())
	args.update(st.init_paths(args))
	return args

def snp_sharing(args, samples):
	""" Quantify snp sharing between sample pairs """
	outfile = open(args['out'], 'w')
	header = ['sample1', 'sample2', 'depth1', 'depth2', 'snps1', 'snps2', 'union', 'xsect']
	outfile.write('\t'.join(header)+'\n')
	counts = {}
	for s1 in samples.values():
		for s2 in samples.values():
			if s1.id == s2.id : continue
			union = len(s1.snps.union(s2.snps))
			xsect = len(s1.snps.intersection(s2.snps))
			rec = [s1.id, s2.id, s1.mean_depth, s2.mean_depth, len(s1.snps), len(s2.snps), union, xsect]
			outfile.write('\t'.join([str(_) for _ in rec])+'\n')

if __name__ == '__main__':
	args = parse_arguments()
	utility.print_copyright()
	
	print("Loading samples...")
	samples = st.init_samples(args)
		
	print("Identifying SNPs...")
	st.identify_snps(args, samples)
	
	print("Counting shared SNPs...")
	snp_sharing(args, samples)
	







