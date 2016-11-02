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
		epilog="""Examples:
1) Run with defaults:
gene_sharing.py --indir /path/to/species --out distances.txt

2) Run a quick test:
gene_sharing.py --indir /path/to/species --out distances.txt --max_genes 1000 --max_samples 10

3) Use a different distance metric:
gene_sharing.py --indir /path/to/species --out distances.txt --distance manhattan

4) Use a lenient cutoff for determining gene presence absence:
gene_sharing.py --indir /path/to/species --out distances.txt --cutoff 0.10

5) Use a strict cutoff for determining gene presence absence:
gene_sharing.py --indir /path/to/species --out distances.txt --cutoff 0.75
""")
	parser.add_argument('--indir', metavar='PATH', type=str, required=True,
		help="""Path to output from `merge_midas.py genes` for one species
Directory should be named according to a species_id and contains files 'genes_*.txt')""")
	parser.add_argument('--out', metavar='PATH', type=str, required=True,
		help="""Path to output file""")
	parser.add_argument('--max_genes', metavar='INT', type=int,
		help="""Maximum number of genes to use. Useful for quick tests (use all)""")
	parser.add_argument('--max_samples', metavar='INT', type=int,
		help="""Maximum number of samples to use. Useful for quick tests (use all)""")
	parser.add_argument('--distance', choices=['jaccard', 'euclidean', 'manhattan'], default='jaccard',
		help="""Metric to use for computing distances (jaccard)""")
	parser.add_argument('--dtype', choices=['presabs', 'copynum'], default='presabs',
		help="""Data type to use for comparing genes (presabs)""")
	parser.add_argument('--cutoff', metavar='FLOAT', type=float, default=0.35,
		help="""Cutoff to use for determining presence absence (0.35)""")

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
			sys.exit("\nError: Input file does not exist: %s\n" % inpath)
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
	
	print ("Reading gene copy-number matrix\n")
	usecols = range(args['max_samples']+1) if args['max_samples'] else None
	data = pd.read_table(args['copynum'], index_col='gene_id', nrows=args['max_genes'], usecols=usecols )
	
	if args['dtype'] == 'presabs':
		print ("Converting to gene presence-absence matrix\n")
		data = data.applymap(cast_presabs)
	
	print ("Computing distances between sample pairs\n")
	# setup output file
	outfile = open(args['out'], 'w')
	fields = ['sample1', 'sample2', 'count1', 'count2', 'count_both', 'count_either', 'distance']
	outfile.write('\t'.join(fields)+'\n')

	# loop over sample pairs
	pairs = itertools.combinations(data.columns.values, 2)
	for index, pair in enumerate(pairs):
		
		sample1, sample2 = pair
		# compute stats
		
		if args['dtype'] == 'presabs':
			set1 = set(data.index[data[sample1]==1]); count1 = len(set1)
			set2 = set(data.index[data[sample2]==1]); count2 = len(set2)
			set_both = set1.intersection(set2); count_both = len(set_both)
			set_either = set1.union(set2); count_either = len(set_either)
		else:
			count1 = sum(data[sample2])
			count2 = sum(data[sample2])
			count_both = sum([min(i,j) for i,j in zip(data[sample1], data[sample2])])
			count_either = sum([max(i,j) for i,j in zip(data[sample1], data[sample2])])

		if args['distance'] == 'jaccard':
			distance = 1-(float(count_both)/count_either) if count_either > 0 else 0
		elif args['distance'] == 'euclidean':
			print 1
			distance = compute_euclidian(data, sample1, sample2)
		else:
			print 2
			distance = compute_manhattan(data, sample1, sample2)
		
		# write stats
		values = [sample1, sample2, count1, count2, count_both, count_either, distance]
		outfile.write('\t'.join([str(_) for _ in values])+'\n')




