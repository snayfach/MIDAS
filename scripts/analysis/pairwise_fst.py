#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip, numpy as np, random, analyze_snps as st
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import utility

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

def pi(x):
	""" Nucleotide diversity """
	return(2*x*(1-x))

def pair_fst(f1, f2):
	f = (f1 + f2)/2
	p = f * (1-f)
	if p == 0:
		return 0, 0
	else:
		p1 = f1 * (1-f1)
		p2 = f2 * (1-f2)
		fst = (p - 0.5 * (p1 + p2) )/p
		return fst, p


if __name__ == '__main__':
	args = parse_arguments()
	utility.print_copyright()
	
	print("Loading samples...")
	samples = st.init_samples(args)
		
	print("Computing Fst...")
	pairs = {}
	for id1 in samples:
		for id2 in samples:
			pairs[id1, id2] = {'fst':0, 'var':0, 'sites':0}

	for site in st.parse_sites(args):

		mean_freq = np.mean([float(freq) for freq in site.freqs if freq != 'NA'])
		if not st.is_snp(mean_freq, args['maf']):
			continue

		for id1, freq1 in zip(site.sample_ids, site.freqs):
			for id2, freq2 in zip(site.sample_ids, site.freqs):
				if freq1 == 'NA' or freq2 == 'NA':
					continue
				else:
					fst, var = pair_fst(float(freq1), float(freq2))
					pairs[id1, id2]['fst'] += fst * var
					pairs[id1, id2]['var'] += var
					pairs[id1, id2]['sites'] += 1

	print("Writing results...")
	outfile = open(args['out'], 'w')
	fields = ['sample1', 'sample2', 'sites', 'fst']
	outfile.write('\t'.join(fields)+'\n')
	for sample1, sample2 in pairs:
		values = pairs[sample1, sample2]
		fst = values['fst']/values['var'] if values['var'] > 0 else 0
		rec = [sample1, sample2, str(values['sites']), str(fst)]
		outfile.write('\t'.join(rec)+'\n')




