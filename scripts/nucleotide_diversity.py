#!/usr/bin/env python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.2'

import argparse, sys, os, gzip, numpy as np

class Sample:
	""" Base class for samples """
	def __init__(self, sample_id):
		self.id = sample_id
		self.total_diversity = 0
		self.count_sites = 0
		self.ref_freq = []

	def compute_pi(self):
		if self.count_sites > 0:
			self.pi = self.total_diversity/self.count_sites
		else:
			self.pi = 'NA'

def print_copyright():
	print ("")
	print ("PhyloCNV: species abundance and strain-level genomic variation from metagenomes")
	print ("version %s; github.com/snayfach/PhyloCNV" % __version__)
	print ("Copyright (C) 2015 Stephen Nayfach")
	print ("Freely distributed under the GNU General Public License (GPLv3)")
	print ("")

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__),
		description="""Estimate core-genome nucleotide diversity of a species within and between samples""")
	parser.add_argument('-i', type=str, dest='input', required=True,
		help="""input to reference frequency matrix""")
	parser.add_argument('-o', dest='output', type=str, required=True,
		help="""path to output file""")
	parser.add_argument('--max_sites', dest='max_sites', type=int, required=False,
		help="""maximum number of sites to use from input matrix. useful for quick tests (use all)""")
	args = vars(parser.parse_args())
	check_args(args)
	return args

def check_args(args):
	""" Check validity of command line arguments """
	pass

def parse_matrix(inpath):
	""" Parse SNP matrix """
	ext = inpath.split('.')[-1]
	infile = gzip.open(inpath) if ext == 'gz' else open(inpath)
	samples = next(infile).rstrip().split()[2:]
	for line in infile:
		values = line.rstrip().split()
		snp = values[0:2]
		freqs = values[2:]
		yield snp, samples, freqs

def fetch_samples(inpath):
	""" Get list of samples from SNP matrix """
	ext = inpath.split('.')[-1]
	infile = gzip.open(inpath) if ext == 'gz' else open(inpath)
	samples = next(infile).rstrip().split()[2:]
	return(samples)

def pi(x):
	""" Nucleotide diversity """
	return(2*x*(1-x))

def init_samples(args):
	""" Initialize samples """
	samples = {}
	for sample_id in fetch_samples(args['input']):
		samples[sample_id] = Sample(sample_id)
	samples['between'] = Sample('between')
	samples['between_cons'] = Sample('between_cons')
	return samples

def compute_nucleotide_diversity(args, samples):
	""" Compute per-bp nucleotide diversity for samples """
	# total diversity
	for index, values in enumerate(parse_matrix(args['input'])):
		f = []
		snp, sample_ids, freqs = values
		if args['max_sites'] and index >= args['max_sites']: break
		for sample_id, freq in zip(sample_ids, freqs):
			if freq == 'NA': continue
			samples[sample_id].total_diversity += pi(float(freq))
			samples[sample_id].count_sites += 1
			f.append(float(freq))
		if len(f) > 1:
			samples['between'].total_diversity += pi(np.mean(f))
			samples['between'].count_sites += 1
			samples['between_cons'].total_diversity += pi(np.mean([round(i) for i in f]))
			samples['between_cons'].count_sites += 1
	# per-bp diversity
	for sample_id in sorted(samples.keys()):
		samples[sample_id].compute_pi()

def write_results(args, samples):
	ext = args['output'].split('.')[-1]
	outfile = gzip.open(args['output'], 'w') if ext == 'gz' else open(args['output'], 'w')
	outfile.write('\t'.join(['sample_id', 'nucleotide_diversity', 'count_sites'])+'\n')
	for sample_id in sorted(samples.keys()):
		sample = samples[sample_id]
		record = [sample.id, sample.pi, sample.count_sites]
		outfile.write('\t'.join([str(x) for x in record])+'\n')

if __name__ == '__main__':
	args = parse_arguments()
	print_copyright()
	print("Estimating nucleotide diversity within and between samples...\n")
	samples = init_samples(args)
	compute_nucleotide_diversity(args, samples)
	write_results(args, samples)

