#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip, numpy as np, pandas as pd, itertools
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import utility

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: gene_sharing.py [options]
""")
	parser.add_argument('--out', metavar='PATH', type=str,
		help="""path to output file""")
	parser.add_argument('--inbase', metavar='PATH', type=str, required=True,
		help="""basename for input files""")
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
	for ext in ['presabs', 'gene_depth', 'copynum']:
		#inpath = '%s.genes.%s' % (args['inbase'], ext)
		inpath = '%s.%s' % (args['inbase'], ext) # TESTING
		if os.path.isfile(inpath):
			paths[ext] = inpath
		else:
			sys.exit('Input file does not exist: %s' % inpath)
	return paths

def euclidian(df, s1, s2):
	return(np.sqrt(sum((df[s1]-df[s2])**2)))

def jaccard(df, s1, s2):
	a = set(df.index[df[s1]==1])
	b = set(df.index[df[s2]==1])
	u = a.union(b)
	x = a.intersection(b)
	return 1 - len(x)/float(len(u))

def manhattan(df, s1, s2):
	return(float(sum(abs(df[s1]-df[s2]))))

def count_diff(df, s1, s2, lower, upper):
	x1 = [0 if _ < lower else 3 if _ > upper else 1 for _ in df[s1]]
	x2 = [0 if _ < lower else 3 if _ > upper else 1 for _ in df[s2]]
	return sum([1 if sum(_) == 3 else 0 for _ in zip(x1, x2)])

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
	outfile = open(args['out'], 'w')
	fields = ['sample1', 'sample2', 'genes1', 'copynum1', 'genes2', 'copynum2', 'count_diff', 'jaccard', 'euclidian', 'manhattan']
	outfile.write('\t'.join(fields)+'\n')
	pairs = itertools.combinations(presabs.columns.values, 2)
	for index, pair in enumerate(pairs):
		id1, id2 = pair
		rec = [id1, id2, sum(presabs[id1]), sum(copynum[id1]), sum(presabs[id2]), sum(copynum[id2])]
		rec.append(count_diff(copynum, id1, id2, args['lower'], args['upper']))
		rec.append(jaccard(presabs, id1, id2))
		rec.append(euclidian(copynum, id1, id2))
		rec.append(manhattan(copynum, id1, id2))
		outfile.write('\t'.join([str(_) for _ in rec])+'\n')

