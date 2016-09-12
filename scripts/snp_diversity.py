#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip, numpy as np, random
from midas import utility
from midas.analyze import snp_matrix

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: 
  estimate strain-level genomic diversity of an individual bacterial species
  diversity estimated across pooled metagenomic samples, or for individual samples
  diversity estimated genome-wide, for individual genes, or for certain types of genomic sites
  
Usage: snp_diversity.py indir [options]

""",
		epilog="""Examples:

""")
	parser.add_argument('--indir', metavar='PATH', type=str, required=True,
		help="""path to output from `merge_midas.py snps` for one species
directory should be named according to a species_id and contains files 'snps_*.txt')""")
	parser.add_argument('--out', metavar='PATH', type=str, required=True,
		help="""path to output file containing list of markers""")

	diversity = parser.add_argument_group("Diversity options")
	diversity.add_argument('--genomic_type', choices=['genome-wide'], default='genome-wide',
		help="""compute diversity for individual genes or genome-wide (genome-wide)""")
	diversity.add_argument('--sample_type', choices=['per-sample', 'pooled-samples'], default='per-sample',
		help="""compute diversity for individual samples or for pooled reads across samples (per-sample)""")
	diversity.add_argument('--site_type', choices=['ALL','NC','CDS','1D','2D','3D','4D'], default='ALL',
		help="""compute diversity using subset of genomic sites sites (ALL)
ALL=all-sites, NC=non-coding, CDS=coding, XD=X-fold-degenerate-sites""")
	diversity.add_argument('--weight_by_depth', action="store_true", default=False,
		help="""weight samples by sequencing depth when --sample_type=pooled-samples""")
	diversity.add_argument('--rand_reads', type=int, metavar='INT',
		help="""randomly select N reads from each sample for each genomic site """)
	diversity.add_argument('--rand_samples', type=int, metavar='INT',
		help="""randomly select N samples from each genomic site""")
	diversity.add_argument('--rand_sites', type=int, metavar='INT',
		help="""randomly select N genomic sites""")
		
	sample = parser.add_argument_group("Sample filters (select subset of samples from INDIR)")
	sample.add_argument('--samples', type=str, metavar='STR', required=False,
		help="""comma-separated list of samples to use for computing diversity metrics""")
	sample.add_argument('--sample_depth', dest='sample_depth', type=float, default=0.0, metavar='FLOAT',
		help="""minimum average read depth per sample (0.0)""")
	sample.add_argument('--fract_cov', dest='fract_cov', type=float, default=0.0, metavar='FLOAT',
		help="""fraction of reference sites covered by at least 1 read (0.0)""")
	sample.add_argument('--max_samples', type=int, metavar='INT', default=float('Inf'),
		help="""maximum number of samples to process.
useful for quick tests (use all)""")

	snps = parser.add_argument_group("Site filters (select subset of genomic sites from INDIR)")
	snps.add_argument('--site_depth', type=int, default=3, metavar='INT',
		help="""minimum number of mapped reads per site (3)
a high value like 20 will result in accurate allele frequencies, but may discard many sites.
a low value like 1 will retain many sites but may not result in accurate allele frequencies""")
	snps.add_argument('--site_prev', type=float, default=0.95, metavar='FLOAT',
		help="""site has at least <site_depth> coverage in at least <site_prev> proportion of samples (0.95)
a value of 1.0 will select sites that have sufficent coverage in all samples.
a value of 0.0 will select all sites, including those with low coverage in many samples 
NAs recorded for included sites with less than <site_depth> in a sample """)
	snps.add_argument('--site_maf', type=float, default=0.0, metavar='FLOAT',
		help="""minimum average-minor-allele-frequency of site across samples (0.0)
setting this above zero (e.g. 0.01, 0.02, 0.05) will only retain variable sites
by default invariant sites are also retained.""")
	snps.add_argument('--site_ratio', type=float, default=float('Inf'), metavar='INT',
		help="""maximum ratio of site-depth to mean-genome-depth (None)
a value of 10 will filter genomic sites with 10x high coverage than the genomic background""")
	snps.add_argument('--max_sites', type=int, default=float('Inf'), metavar='INT',
		help="""maximum number of sites to include in output (use all)
useful for quick tests""")
	args = vars(parser.parse_args())
	if args['site_type'] == 'ALL': # format site_type
		args['site_type'] = ['NC','1D','2D','3D','4D']
	elif args['site_type'] == 'CDS':
		args['site_type'] = ['1D','2D','3D','4D']
	else:
		args['site_type'] = [args['site_type']]
	return args

def print_args(args):
	lines = []
	lines.append("===========Parameters===========")
	lines.append("Script: snp_diversity.py")
	lines.append("Input directory: %s" % args['indir'])
	lines.append("Output file: %s" % args['out'])
	lines.append("Diversity options:")
	lines.append("   -genomic_type: %s" % args['genomic_type'])
	lines.append("   -sample_type: %s" % args['sample_type'])
	lines.append("   -site_type: %s" % args['site_type'])
	lines.append("   -weight_by_depth: %s" % args['weight_by_depth'])
	lines.append("   -rand_reads: %s" % args['rand_reads'])
	lines.append("   -rand_sites: %s" % args['rand_sites'])
	lines.append("Sample filters:")
	lines.append("   --samples: %s" % args['samples'])
	lines.append("   --sample_depth: %s" % args['sample_depth'])
	lines.append("   --fract_cov: %s" % args['fract_cov'])
	lines.append("   --max_samples: %s" % args['max_samples'])
	lines.append("Site filters:")
	lines.append("   --site_depth: %s" % args['site_depth'])
	lines.append("   --site_prev: %s" % args['site_prev'])
	lines.append("   --site_maf: %s" % args['site_maf'])
	lines.append("   --site_ratio: %s" % args['site_ratio'])
	lines.append("   --max_sites: %s" % args['site_ratio'])
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
	if args['rand_reads'] > args['site_depth']:
		sys.exit("\nError: --rand_reads cannot exceed --site_depth")

class Diversity:
	def __init__(self):
		self.sites = 0
		self.samples = 0
		self.snps = 0
		self.pi = 0

def init_pi(args, samples):
	pi = {}
	if args['sample_type'] == 'per-sample':
		for s in samples:
			if s.pass_qc:
				if args['genomic_type'] == 'genome-wide':
					pi[s.id] = Diversity()
				else:
					pi[s.id] = {}
					for g in genes:
						pi[s.id][g.id] = Diversity()
	else:
		if args['genomic_type'] == 'genome-wide':
			pi = Diversity()
		else:
			for g in genes:
				pi[s.id][g.id] = Diversity()
	return pi

def maf(x):
	return x if x <= 0.5 else 1-x

def compute_pi(args, samples):
	pi = init_pi(args, samples)
	index = 0
	for site in snp_matrix.parse_sites(args['indir'], samples): # loop over sites
		
		# stop early
		if index >= args['max_sites']:
			break
			
		# prune low quality samples
		site.prune_samples(args['site_depth'], args['site_ratio'])
		
		# compute site summary stats
		site.summary_stats(args['weight_by_depth'])
		
		# filter site
		if site.filter_site(args['site_prev'], args['site_maf'], args['site_type']):
			continue
		else:
			index += 1

		# downsample samples
		if args['rand_samples']:
			site.resample_samples(args['rand_samples'])
		
		# downsample reads
		if args['rand_reads'] and site.maf > 0:
			site.resample_reads(args['rand_reads'])
		
		# compute pi for pooled-samples
		if args['sample_type'] == 'pooled-samples':
			site_pi = site.pooled_pi(args['site_depth'], args['weight_by_depth'])
			pi.pi += site_pi
			pi.sites += 1

		# compute pi for per-samples
		else:
			for sample, values in site.sample_values().items():
				if values['depth'] < args['site_depth']:
					continue
				if args['rand_reads'] and maf(values['ref_freq']) > 0: # downsample reads
					count_ref = int(round(values['ref_freq'] * values['depth']))
					count_alt = values['depth'] - count_ref
					alleles = random.sample([1] * count_ref + [0] * count_alt, args['rand_reads'])
					values['ref_freq'] = np.mean(alleles)
				site_maf = maf(values['ref_freq'])
				pi[sample.id].pi += 2 * site_maf * (1-site_maf)
				pi[sample.id].sites += 1
	return pi

def write_pi(args, samples, pi):
	outfile = open(args['out'], 'w')
	if args['sample_type'] == 'pooled-samples':
		h = ['count_sites', 'diversity']
		outfile.write('\t'.join([str(_) for _ in h])+'\n')
		r = [pi.sites, pi.pi]
		outfile.write('\t'.join([str(_) for _ in r])+'\n')
	else:
		h = ['sample_id', 'count_sites', 'diversity']
		outfile.write('\t'.join([str(_) for _ in h])+'\n')
		for s in samples:
			if s.pass_qc:
				r = [s.id, pi[s.id].sites, pi[s.id].pi]
				outfile.write('\t'.join([str(_) for _ in r])+'\n')

if __name__ == '__main__':
	args = parse_arguments()
	check_args(args)
	utility.print_copyright()
	print_args(args)
	
	print("Selecting subset of samples...")
	samples = snp_matrix.fetch_samples(args)
	print("   %s samples selected") % len([_ for _ in samples if _.pass_qc])
	
	print("Estimating diversity metrics...")
	pi = compute_pi(args, samples)

	print("Writing results to output file...")
	write_pi(args, samples, pi)


