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
Usage: count_snps.py outdir [options]
""")
	parser.add_argument('outdir', type=str,
		help="""path to output directory""")
	parser.add_argument('--ref_freq', metavar='PATH', type=str, required=True,
		help="""input to allele frequency matrix""")
	parser.add_argument('--max_sites', metavar='INT', type=int, required=False,
		help="""maximum number of sites to use from input matrix
useful for quick tests (use all)""")
	parser.add_argument('--maf', type=float, metavar='FLOAT', required=False, default=0.05,
		help="""minor allele frequency for SNP definition (0.05)""")
	parser.add_argument('--iter', type=int, metavar='INT', required=False, default=10,
		help="""number of rarefactions fo SNP counting (10)""")
	parser.add_argument('--samples', type=str, metavar='STR', required=True,
		help="""comma-separated list of sample ids (none)""")

	args = vars(parser.parse_args())
	return args

def parse_matrix(inpath):
	""" Parse SNP matrix """
	ext = inpath.split('.')[-1]
	infile = gzip.open(inpath) if ext == 'gz' else open(inpath)
	samples = next(infile).rstrip().split()[2:]
	for line in infile:
		rec = line.rstrip().split()
		snp = rec[0:2]
		values = rec[2:]
		yield snp, samples, values

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
	sample_ids = args['samples'].split(',')
	for sample_id in fetch_samples(args['ref_freq']):
		if sample_id in sample_ids:
			samples[sample_id] = Sample(sample_id)
	samples['between'] = Sample('between')
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
	for snp, sample_ids, freqs in parse_matrix(args['ref_freq']):
		for sample_id, freq in zip(sample_ids, freqs):
			if sample_id not in samples:
				continue
			elif freq == 'NA':
				continue
			elif is_snp(float(freq), args['maf']):
				samples[sample_id].snps.add('|'.join(snp))
			samples[sample_id].sites += 1
		if len(freqs) > 1:
			freq = np.mean([float(freq) for freq in freqs if freq != 'NA'])
		if args['max_sites'] and index >= args['max_sites']:
			break
		else:
			index += 1

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
	outfile = open('%s/rarefaction.txt' % args['outdir'], 'w')
	outfile.write('count_samples\titer\tcount_union\tcount_xsect\n')
	for n in counts:
		union = counts[n]['union']
		xsect = counts[n]['xsect']
		for i, j, k in zip(range(args['iter']), union, xsect):
			r = [n, i, j, k]
			outfile.write('\t'.join([str(x) for x in r])+'\n')


def compute_prevalence(args, samples):
	snps = set([])
	for sample in samples.values():
		snps |= sample.snps
	snps = dict([(_,0) for _ in snps])
	for sample in samples.values():
		for snp in sample.snps:
			snps[snp] += 1
	return snps

def write_prevalence(args, prevalence, samples):
	outfile = open('%s/prevalence.txt' % args['outdir'], 'w')
	outfile.write('snp_id\tcount_samples\tfraction_samples\n')
	counts = prevalence.items()
	counts.sort(key=lambda x: x[1], reverse=True)
	for snp_id, count in counts:
		fraction = float(count)/len(samples)
		r = [snp_id, count, fraction]
		outfile.write('\t'.join([str(x) for x in r])+'\n')

if __name__ == '__main__':
	args = parse_arguments()
	
	print("Loading samples...")
	samples = init_samples(args)

	print("Identifying SNPs...")
	identify_snps(args, samples)
	#write_snps(args, counts)
	
	print("Computing SNP prevalence...")
	prevalence = compute_prevalence(args, samples)
	write_prevalence(args, prevalence, samples)
		
	print("Performing SNP rarefaction...")
	counts = count_snps(args, samples)
	write_counts(args, counts)

	#print("Computing pairwise statistics...")








