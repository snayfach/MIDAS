#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip, numpy as np, pandas as pd, itertools
from midas import utility

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""Description:
This script will compare the gene content between all pairs of samples
Gene presence-absence and gene copy number are used for comparison

Usage: gene_sharing.py indir outfile [options]
""",
		epilog="""Output format:
1) sample1: first sample identifier
2) sample2: second sample identifier
3) genes1: number of present genes in sample1 (see --cutoff)
5) genes2: number of present genes in sample2 (see --cutoff)
5) genes_union: number of present genes in sample2 (see --cutoff)
5) genes_shared: number of present genes in sample2 (see --cutoff)
8) jaccard: dissimilarity between gene sets (see --cutoff)

4) copynum1: total gene copy number in sample1
6) copynum2: total gene copy number in sample2
9) euclidian: euclidian distance between gene copy numbers
10) manhattan: manhattan distance between gene copy numbers
""")
	parser.add_argument('indir', metavar='PATH', type=str,
		help="""directory for input files""")
	parser.add_argument('outfile', metavar='PATH', type=str,
		help="""path to output file""")
	parser.add_argument('--cutoff', metavar='FLOAT', type=float, default=0.35,
		help="""cutoff to use for determining\ngene presence-absence (0.35)""")
	parser.add_argument('--lower', metavar='FLOAT', type=float, default=0.05,
		help="""lower cutoff for confident classification of gene absence (0.05)""")
	parser.add_argument('--upper', metavar='FLOAT', type=float, default=0.50,
		help="""upper cutoff for confident classification of gene presence (0.50)""")
	parser.add_argument('--max_genes', metavar='INT', type=int,
		help="""maximum number of genes to use. useful for quick tests (use all)""")
	parser.add_argument('--max_samples', metavar='INT', type=int,
		help="""maximum number of samples to use. useful for quick tests (use all)""")
	args = vars(parser.parse_args())
	args.update(init_paths(args))
	return args

def init_paths(args):
	paths = {}
	for ext in ['presabs', 'depth', 'copynum']:
		inpath = '%s/genes_%s.txt' % (args['indir'], ext)
		if os.path.isfile(inpath):
			paths[ext] = inpath
		else:
			sys.exit('Input file does not exist: %s' % inpath)
	return paths

def compute_euclidian(df, s1, s2):
	return(np.sqrt(sum((df[s1]-df[s2])**2)))

def compute_jaccard(df, s1, s2, type):
	if type=='binary':
		a = set(df.index[df[s1]==1])
		b = set(df.index[df[s2]==1])
		x = a.intersection(b)
		u = a.union(b)
		j = 1-(len(x)/float(len(u))) if len(u) > 0 else 'NA'
		return (len(a), len(b), len(x), len(u), j)
	else:
		a = sum(df[s1])
		b = sum(df[s2])
		x = sum([min(i,j) for i,j in zip(df[s1], df[s2])])
		u = sum([max(i,j) for i,j in zip(df[s1], df[s2])])
		j = 1-(x/u) if u > 0 else 'NA'
		return (a, b, x, u, j)

def compute_manhattan(df, s1, s2):
	return(float(sum(abs(df[s1]-df[s2]))))

def count_fixed_diffs(df, s1, s2, lower, upper):
	x1 = [0 if _ < lower else 3 if _ > upper else 1 for _ in df[s1]]
	x2 = [0 if _ < lower else 3 if _ > upper else 1 for _ in df[s2]]
	a = sum([1 for _ in df[s1] if _ > upper ]) # a present
	b = sum([1 for _ in df[s2] if _ > upper ]) # b present
	a_not_b = sum([1 for i,j in zip(df[s1], df[s2]) if i>upper and j<lower]) # a but not b present
	b_not_a = sum([1 for i,j in zip(df[s1], df[s2]) if i<lower and j>upper]) # b but not a present
	return a, b, a_not_b, b_not_a

def cast_presabs(x):
	""" Applies cutoff to x; args must be global """
	return 1 if x > args['cutoff'] else 0

if __name__ == '__main__':
	
	args = parse_arguments()
	
	print ("Reading gene copy-number matrix...")
	usecols = range(args['max_samples']+1) if args['max_samples'] else None
	copynum = pd.read_table(args['copynum'], index_col='gene_id', nrows=args['max_genes'], usecols=usecols )
	print("  %s Gb maximum memory") % utility.max_mem_usage()
	
	print ("Generating gene presence-absence matrix ...")
	presabs = copynum.applymap(cast_presabs)
	print("  %s Gb maximum memory") % utility.max_mem_usage()
	
	print ("Computing distances between sample pairs ...")
	outfile = open(args['outfile'], 'w')
	fields = ['sample1', 'sample2',
			  'upper1', 'upper2', 'upper1_lower2', 'upper2_lower1',
	          'present1', 'present2', 'present_both', 'present_either', 'jaccard_presabs',
			  'copynum1', 'copynum2', 'copynum_both', 'copynum_either', 'jaccard_copynum',
			  'euclidian_copynum', 'manhattan_copynum']
	outfile.write('\t'.join(fields)+'\n')
	pairs = itertools.combinations(presabs.columns.values, 2)
	
	for index, pair in enumerate(pairs):
		# compute stats
		sample1, sample2 = pair
		upper1, upper2, upper1_lower2, upper2_lower1 = count_fixed_diffs(copynum, sample1, sample2, args['lower'], args['upper'])
		present1, present2, present_both, present_either, jaccard_presabs = compute_jaccard(presabs, sample1, sample2, 'binary')
		copynum1, copynum2, copynum_both, copynum_either, jaccard_copynum = compute_jaccard(copynum, sample1, sample2, 'presabs')
		euclidian_copynum = compute_euclidian(copynum, sample1, sample2)
		manhattan_copynum = compute_manhattan(copynum, sample1, sample2)
		
		# write stats
		values = [sample1, sample2,
				  upper1, upper2, upper1_lower2, upper2_lower1,
				  present1, present2, present_both, present_either, jaccard_presabs,
				  copynum1, copynum2, copynum_both, copynum_either, jaccard_copynum,
				  euclidian_copynum, manhattan_copynum]
		outfile.write('\t'.join([str(_) for _ in values])+'\n')




