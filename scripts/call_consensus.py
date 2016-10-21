#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip, numpy as np, random, csv
from midas.utility import print_copyright
from midas import utility, parse

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description:
Build FASTA file of consensus sequences for a species per sample
Useful for building phylogenetic trees
Before running this script, you'll need to have run `merge_midas.py snps`

Usage: call_consensus.py --indir <PATH> --out <PATH> [options]
""",
		epilog="""
Examples:
1) Build multi-FASTA of core-genome sequences (recommended)
-core-genome sites defined as >=5 reads in >=90% of samples
-use only variable positions (>=1% minor allele frequency across samples)
-only include samples with sufficient data (>=10x mean-depth, >=40% of sites with >=1 mapped read)
-exclude sites with abnormal depth (>5x mean-depth or <1/5 mean-depth)

call_consensus.py --indir /path/to/snps --out /path/to/seqs --site_maf 0.01 --site_depth 5 --site_prev 0.90 --sample_depth 10 --sample_cov 0.40 --site_ratio 5.0

2) Build multi-FASTA using defaults
call_consensus.py --indir /path/to/snps --out /path/to/seqs

3) Run a quick test
call_consensus.py --indir /path/to/snps --out /path/to/seqs --max_sites 10000

""")
	parser.add_argument('--indir', metavar='PATH', type=str, required=True,
		help="""path to output from `merge_midas.py snps` for one species
directory should be named according to a species_id and contains files 'snps_*.txt')""")
	parser.add_argument('--out', metavar='PATH', type=str, required=True,
		help="""path to output file""")

	sample = parser.add_argument_group("Sample filters (select subset of samples from INDIR)")
	sample.add_argument('--sample_depth', dest='sample_depth', type=float, default=0.0, metavar='FLOAT',
		help="""minimum average read depth per sample (0.0)""")
	sample.add_argument('--sample_cov', dest='fract_cov', type=float, default=0.0, metavar='FLOAT',
		help="""fraction of reference sites covered by at least 1 read (0.0)""")
	sample.add_argument('--max_samples', type=int, metavar='INT', default=float('Inf'),
		help="""maximum number of samples to process.
useful for quick tests (use all)""")
	sample.add_argument('--keep_samples', type=str, metavar='STR', required=False,
		help="""comma-separated list of samples to include
samples will still be subject to other filters""")
	sample.add_argument('--exclude_samples', type=str, metavar='STR', required=False,
		help="""comma-separated list of samples to exclude.
samples will still be subject to other filters""")

	snps = parser.add_argument_group("Site filters (select subset of genomic sites from INDIR)")
	snps.add_argument('--site_depth', type=int, default=2, metavar='INT',
		help="""minimum number of mapped reads per site (2)""")
	snps.add_argument('--site_prev', type=float, default=0.0, metavar='FLOAT',
		help="""site has at least <site_depth> coverage in at least <site_prev> proportion of samples (0.0)
a value of 1.0 will select sites that have sufficent coverage in all samples.
a value of 0.0 will select all sites, including those with low coverage in many samples 
NAs recorded for included sites with less than <site_depth> in a sample """)
	snps.add_argument('--site_maf', type=float, default=0.0, metavar='FLOAT',
		help="""minimum average-minor-allele-frequency of site across samples (0.0)
setting this above zero (e.g. 0.01, 0.02, 0.05) will only retain variable sites
by default invariant sites are also retained.""")
	snps.add_argument('--site_ratio', type=float, default=float('Inf'), metavar='FLOAT',
		help="""maximum ratio of site-depth to mean-genome-depth (None)
a value of 10 will filter genomic sites with 10x high coverage than the genomic background""")
	snps.add_argument('--max_sites', type=int, default=float('Inf'), metavar='INT',
		help="""maximum number of sites to include in output (use all)
useful for quick tests""")
	args = vars(parser.parse_args())
	keep = args['keep_samples'].rstrip(',').split(',') if args['keep_samples'] else None
	exclude = args['exclude_samples'].rstrip(',').split(',') if args['exclude_samples'] else None
	return args

def print_args(args):
	lines = []
	lines.append("===========Parameters===========")
	lines.append("Script: call_consensus.py")
	lines.append("Input directory: %s" % args['indir'])
	lines.append("Output file: %s" % args['out'])
	lines.append("Sample filters:")
	lines.append("  sample_depth: %s" % args['sample_depth'])
	lines.append("  fract_cov: %s" % args['fract_cov'])
	lines.append("  max_samples: %s" % args['max_samples'])
	lines.append("  keep_samples: %s" % args['keep_samples'])
	lines.append("  exclude_samples: %s" % args['exclude_samples'])
	lines.append("Site filters:")
	lines.append("  site_depth: %s" % args['site_depth'])
	lines.append("  site_prev: %s" % args['site_prev'])
	lines.append("  site_maf: %s" % args['site_maf'])
	lines.append("  site_ratio: %s" % args['site_ratio'])
	lines.append("  max_sites: %s" % args['max_sites'])
	sys.stdout.write('\n'.join(lines)+'\n')

def check_args(args):
	if not os.path.isdir(args['indir']):
		sys.exit("Specified input directory '%s' does not exist" % args['indir'])
	if args['site_depth'] < 2:
		sys.exit("\nError: --site_depth must be >=2 to calculate nucleotide variation")
	if args['max_sites'] < 1:
		sys.exit("\nError: --max_sites must be >= 1 to calculate nucleotide variation")
	if args['max_samples'] < 1:
		sys.exit("\nError: --max_samples must be >= 1 to calculate nucleotide variation")
	if args['site_ratio'] < 0:
		sys.exit("\nError: --site_ratio cannot be a negative number")
	if args['site_depth'] < 0:
		sys.exit("\nError: --site_depth cannot be a negative number")
	if args['sample_depth'] < 0:
		sys.exit("\nError: --sample_depth cannot be a negative number")
	if not 0 <= args['site_maf'] <= 1:
		sys.exit("\nError: --site_maf must be between 0 and 1")
	if not 0 <= args['site_prev'] <= 1:
		sys.exit("\nError: --site_prev must be between 0 and 1")
	if not 0 <= args['fract_cov'] <= 1:
		sys.exit("\nError: --fract_cov must be between 0 and 1")

def init_sequences(samples):
	for sample in samples.values():
		sample.seq = ""

def call_consensus(site, sample):
	if not sample.keep:
		return '-' # this could be a parameter
	if sample.depth == 0:
		return '-'
	elif sample.ref_freq >= 0.5:
		return site.ref_allele
	else:
		return sample.alt_allele

def percent_missing(seq):
	return round(100*seq.count('-')/float(len(seq)),2)

def sequence_description(sample):
	desc = {}
	desc['length'] = len(sample.seq)
	desc['percent_missing'] = percent_missing(sample.seq)
	desc['mean_depth'] = round(sample.sample_depth, 2)
	return desc

def write_consensus(args, samples):
	outfile=open(args['out'], 'w')
	for sample_id in sorted(samples):
		sample = samples[sample_id]
		desc = sequence_description(sample)
		outfile.write('>'+sample.id+'\t')
		outfile.write(' '.join(['%s=%s' % (key, value) for key, value in desc.items()])+'\n')
		outfile.write(sample.seq+'\n')
	outfile.close()

if __name__ == '__main__':
	
	# setup args
	args = parse_arguments()
	check_args(args)
	print_copyright()
	print_args(args)
	
	# init species, samples, sequences
	species = parse.Species(args['indir'])
	samples = parse.fetch_samples(
					species,
				    args['sample_depth'],
					args['fract_cov'],
					args['max_samples'],
					args['keep_samples'],
					args['exclude_samples'])
	if len(samples) == 0:
		sys.exit("\nError: no samples satisfied your selection criteria.\nTry running again with more lenient parameters\n")
	
	init_sequences(samples)
	
	# loop over genomic sites
	sites = parse.fetch_sites(species, samples)
	for index, site in enumerate(sites):

		# stop early
		if index >= args['max_sites']: break
			
		# filter low quality samples at site
		site.filter_samples(site_depth=args['site_depth'], site_ratio=args['site_ratio'])

		# compute site summary stats
		site.summary_stats()
		
		# filter site
		site.filter_site(site_prev=args['site_prev'], site_maf=args['site_maf'])
		
		# store consensus
		if site.keep:
			for sample in site.samples.values():
				samples[sample.id].seq += call_consensus(site, sample)

	# write consensus
	write_consensus(args, samples)

