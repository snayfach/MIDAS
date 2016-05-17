#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

# DO:
#
# add distance metric for floating point values in CNV matrices
# add allele sharing metrics
#  -jaccard + min allele frequency
#  -distance metric for floating point allele frequencies

import os, sys

def parse_arguments():
	import argparse
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
	parser.add_argument('-i', dest='in', type=str, help='Input matrix or tree. Accepted file types: .presabs, .copynum, .ref_freq, .tree)', required=True)
	parser.add_argument('-o', dest='out', type=str, help='Output distance matrix', required=True)
	return vars(parser.parse_args())

def jaccard(a,b):
	""" Jaccard dissimilarity between sets a and b """
	union = sum([max(i,j) for i,j in zip(a,b)])
	xsect = sum([min(i,j) for i,j in zip(a,b)])
	return(1 - xsect/float(union))

def read_cnv_matrix(inpath):
	""" Read matrix into memory """
	infile = open(inpath)
	samples = next(infile).rstrip().split()[1:]
	matrix = dict([(s,[]) for s in samples])
	for line in infile:
		values = line.rstrip().split()[1:]
		for sample, value in zip(samples, values):
			matrix[sample].append(int(value))
	return matrix

def jaccard_distmat(matrix):
	""" Compute gene sharing between pairs of samples """
	samples = matrix.keys()
	dist = dict([(s,[]) for s in samples])
	for s1 in samples:
		for s2 in samples:
			a1 = matrix[s1]
			a2 = matrix[s2]
			dist[s1].append(jaccard(a1, a2))
	return (samples, dist)

def cnv_distances(inpath, outpath):
	""" Pairwise distances between samples based on gene presence/absence """
	matrix = read_cnv_matrix(inpath)
	samples, dist = jaccard_distmat(matrix)
	outfile = open(outpath, 'w')
	header = outfile.write('\t%s' % '\t'.join(samples)+'\n')
	for sample in samples:
		outfile.write(sample)
		for value in dist[sample]:
			outfile.write('\t%s' % str(round(value,3)))
		outfile.write('\n')

def phylo_distances(inpath, outpath):
	""" Pairwise distances between samples based on phylogenetic tree """
	import dendropy
	tree = dendropy.Tree.get_from_path(inpath, 'newick')
	dist = dendropy.calculate.treemeasure.PatristicDistanceMatrix(tree)
	samples = tree.taxon_namespace
	outfile = open(outpath, 'w')
	header = outfile.write('\t%s' % '\t'.join([x.label.replace(' ', '_') for x in samples])+'\n')
	for s1 in samples:
		outfile.write(s1.label.replace(' ', '_'))
		for s2 in samples:
			outfile.write('\t%s' % str(round(dist(s1, s2),3)))
		outfile.write('\n')

if __name__ == '__main__':
	
	args = parse_arguments()
	ext = args['in'].split('.')[-1]
	
	if ext == 'presabs':
		cnv_distances(args['in'], args['out'])
	elif ext == 'tree':
		phylo_distances(inpath, outpath)
	else:
		sys.exit("Unrecognized file extension: '%s'\nMust be one of {.presabs, .copynum, .ref_freq, .tree}" % ext)



