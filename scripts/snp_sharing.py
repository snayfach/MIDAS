#!/usr/bin/env python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

# Next: allow user to specify subset of samples in input
#

import argparse, sys, os, gzip, numpy as np, random

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
Usage: snp_sharing.py [options]
""")
	parser.add_argument('--out', metavar='PATH', type=str,
		help="""path to output file""")
	parser.add_argument('--ref_freq', metavar='PATH', type=str, required=True,
		help="""input to allele frequency matrix""")
	parser.add_argument('--depth', metavar='PATH', type=str, required=True,
		help="""input to site depth matrix""")
	parser.add_argument('--summary', metavar='PATH', type=str, required=True,
		help="""input to summary stats""")
	parser.add_argument('--max_sites', metavar='INT', type=int, required=False,
		help="""maximum number of sites to use from input matrix
useful for quick tests (use all)""")
	parser.add_argument('--min_ratio', type=float, metavar='FLOAT', required=False, default=float('-inf'),
		help="""minimum ratio of per-site depth versus genome-average (-inf)""")
	parser.add_argument('--max_ratio', type=float, metavar='FLOAT', required=False, default=float('inf'),
		help="""maximum ratio of per-site depth versus genome-average (inf)""")
	parser.add_argument('--maf', type=float, metavar='FLOAT', required=False, default=0.05,
		help="""minor allele frequency for SNP definition (0.05)""")

	args = vars(parser.parse_args())
	return args

def parse_matrix(inpath, fun):
	""" Parse SNP matrix """
	ext = inpath.split('.')[-1]
	infile = gzip.open(inpath) if ext == 'gz' else open(inpath)
	samples = next(infile).rstrip().split()[2:]
	for line in infile:
		rec = line.rstrip().split()
		snp = rec[0:2]
		#values = [fun(_) for _ in rec[2:] if _ != 'NA']
		yield snp, samples, rec[2:]

def parse_summary(inpath, fun):
	""" Parse SNP matrix """
	ext = inpath.split('.')[-1]
	infile = gzip.open(inpath) if ext == 'gz' else open(inpath)
	samples = next(infile).rstrip().split()[2:]
	for line in infile:
		rec = line.rstrip().split()
		snp = rec[0:2]
		values = [fun(_) for _ in rec[2:]]
		yield snp, samples, values

def fetch_samples(inpath):
	""" Get list of samples from SNP matrix """
	ext = inpath.split('.')[-1]
	infile = gzip.open(inpath) if ext == 'gz' else open(inpath)
	samples = next(infile).rstrip().split()[2:]
	return(samples)

def init_samples(args):
	""" Initialize samples """
	samples = {}
	for sample_id in fetch_samples(args['ref_freq']):
		samples[sample_id] = Sample(sample_id)
	# update with mean-depth
	infile = open(args['summary'])
	fields = next(infile).rstrip().split()
	for line in infile:
		values = line.rstrip().split()
		rec = dict([(f,v) for f,v in zip(fields, values)])
		if rec['run_accession'] in samples:
			sample = samples[rec['run_accession']]
			sample.mean_depth = float(rec['average_depth'])
			sample.min_depth = max(1, args['min_ratio'] * sample.mean_depth)
			sample.max_depth = args['max_ratio'] * sample.mean_depth
	return samples

def is_snp(freq, min_maf):
	""" Determine if freq classified as a SNP or not """
	if min(freq, 1-freq) >= min_maf:
		return True
	else:
		return False

def identify_snps(args, samples):
	""" Store snps per sample """
	index = 0
	ref_freq = parse_matrix(args['ref_freq'], float)
	ref_depth = parse_matrix(args['depth'], int)
	while True:
		try:
			snp, sample_ids, freqs = next(ref_freq)
			snp, sample_ids, counts = next(ref_depth)
		except StopIteration:
			break
		for sample_id, freq, count in zip(sample_ids, freqs, counts):
			# is valid site?
			if (int(count) < samples[sample_id].min_depth
					or int(count) > samples[sample_id].max_depth):
				continue
			else:
				samples[sample_id].sites += 1
			# is SNP?
			if is_snp(float(freq), args['maf']):
				samples[sample_id].snps.add('|'.join(snp))
		if args['max_sites'] and index >= args['max_sites']:
			break
		else:
			index += 1

def snp_sharing(args, samples):
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

	print("Loading samples...\n")
	samples = init_samples(args)
		
	print("Identifying SNPs...\n")
	identify_snps(args, samples)
	
	print("Counting shared SNPs...\n")
	snp_sharing(args, samples)
	







