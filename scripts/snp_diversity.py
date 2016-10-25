#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip, numpy as np, random
from midas.utility import print_copyright, parse_file
from midas.analyze import diversity

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description:
Quantify the genomic diversity of a bacterial population
Diversity computed genome-wide, for different site classes, or for individual genes
Diversity computed for individual metagenomic samples for data pooled across samples
Before running these scripts, you'll need to have run `merge_midas.py snps`

Usage: snp_diversity.py --indir <PATH> --out <PATH> [options]
""",
		epilog="""Examples:

""")
	parser.add_argument('--indir', metavar='PATH', type=str, required=True,
		help="""path to output from `merge_midas.py snps` for one species
directory should be named according to a species_id and contains files 'snps_*.txt')""")
	parser.add_argument('--out', metavar='PATH', type=str, required=True,
		help="""path to output file""")

	diversity = parser.add_argument_group("Diversity options")
	diversity.add_argument('--genomic_type', choices=['genome-wide', 'per-gene'], default='genome-wide',
		help="""compute diversity for individual genes or genome-wide (genome-wide)""")
	diversity.add_argument('--sample_type', choices=['per-sample', 'pooled-samples'], default='per-sample',
		help="""compute diversity for individual samples or for pooled reads across samples (per-sample)""")
	diversity.add_argument('--site_type', choices=['ALL','NC','CDS','1D','2D','3D','4D'], default='ALL',
		help="""compute diversity using subset of genomic sites sites (ALL)
ALL=all-sites, NC=non-coding, CDS=coding, XD=X-fold-degenerate-sites""")
	diversity.add_argument('--weight_by_depth', action="store_true", default=False,
		help="""weight data from samples by sequencing depth when --sample_type=pooled-samples""")
	diversity.add_argument('--rand_reads', type=int, metavar='INT',
		help="""randomly select N reads from each sample for each genomic site """)
	diversity.add_argument('--replace_reads', action='store_true', default=False,
		help="""reads drawn with replacement""")
	diversity.add_argument('--rand_samples', type=int, metavar='INT',
		help="""randomly select N samples from each genomic site""")
	diversity.add_argument('--rand_sites', type=float, metavar='FLOAT',
		help="""randomly select X proportion of high-quality genomic sites""")

	sample = parser.add_argument_group("Sample filters (select subset of samples from INDIR)")
	sample.add_argument('--sample_depth', dest='sample_depth', type=float, default=0.0, metavar='FLOAT',
		help="""minimum average read depth per sample (0.0)""")
	sample.add_argument('--sample_cov', dest='fract_cov', type=float, default=0.0, metavar='FLOAT',
		help="""fraction of reference sites covered by at least 1 read (0.0)""")
	sample.add_argument('--max_samples', type=int, metavar='INT', default=float('Inf'),
		help="""maximum number of samples to process.
useful for quick tests (use all)""")
	sample.add_argument('--keep_samples', type=str, metavar='STR', required=False,
		help="""comma-separated list of samples to use for computing diversity metrics.
samples will still be subject to other filters""")
	sample.add_argument('--exclude_samples', type=str, metavar='STR', required=False,
		help="""comma-separated list of samples to exclude from computing diversity metrics.
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
	snps.add_argument('--site_freq', type=float, default=0.0, metavar='FLOAT',
		help="""minimum combined frequency of reference allele and major alternate allele across samples (0.0)
most mapped reads should match the reference allele or major alternate allele
set this value to exclude sites with multiple alternate alleles
for example, --site_freq=0.95 excludes these sites:
  A(ref)=0.80, T(alt1)=0.10, C(alt2)=0.05, G(alt3)=0.05 (Freq_A + Freq_T = 0.90 < 0.95)
  A(ref)=0.00, C(alt1)=0.80, G(alt2)=0.20, T(alt3)=0.00 (Freq_A + Freq_C = 0.80 < 0.95)""")
	snps.add_argument('--max_sites', type=int, default=float('Inf'), metavar='INT',
		help="""maximum number of sites to include in output (use all)
useful for quick tests""")
	args = vars(parser.parse_args())
	format_site_type(args)
	format_sample_lists(args)
	return args

def format_sample_lists(args):
	keep = args['keep_samples'].rstrip(',').split(',') if args['keep_samples'] else None
	exclude = args['exclude_samples'].rstrip(',').split(',') if args['exclude_samples'] else None

def format_site_type(args):
	if args['site_type'] == 'ALL':
		args['site_type'] = ['NC','1D','2D','3D','4D']
	elif args['site_type'] == 'CDS':
		args['site_type'] = ['1D','2D','3D','4D']
	else:
		args['site_type'] = [args['site_type']]

def print_args(args):
	lines = []
	lines.append("===========Parameters===========")
	lines.append("Command: %s" % ' '.join(sys.argv))
	lines.append("Script: snp_diversity.py")
	lines.append("Input directory: %s" % args['indir'])
	lines.append("Output file: %s" % args['out'])
	lines.append("Diversity options:")
	lines.append("  genomic_type: %s" % args['genomic_type'])
	lines.append("  sample_type: %s" % args['sample_type'])
	lines.append("  site_type: %s" % args['site_type'])
	lines.append("  weight_by_depth: %s" % args['weight_by_depth'])
	lines.append("  rand_reads: %s" % args['rand_reads'])
	lines.append("  replace_reads: %s" % args['replace_reads'])
	lines.append("  rand_samples: %s" % args['rand_samples'])
	lines.append("  rand_sites: %s" % args['rand_sites'])
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
	lines.append("  site_freq: %s" % args['site_freq'])
	lines.append("  max_sites: %s" % args['max_sites'])
	print ("===============================")
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
	if args['rand_reads'] > args['site_depth'] and not args['replace_reads']:
		sys.exit("\nError: --rand_reads cannot exceed --site_depth when --replace_reads=False")
	if args['rand_sites'] and (args['rand_sites'] < 0 or args['rand_sites'] > 1):
		sys.exit("\nError: --rand_sites must be between 0 and 1")
	if 'NC' in args['site_type'] and args['genomic_type'] == 'per-gene':
		sys.exit("\nError: --site_type cannot be NC if --genomic_type is per-gene")
	# set min # of rand reads

class Diversity:
	def __init__(self):
		self.sites = 0
		self.samples = 0
		self.snps = 0
		self.pi = 0
		self.depth = 0

def list_genes(args):
	""" List the set of genes from snps_info.txt """
	genes = set([])
	for r in parse_file('%s/snps_info.txt' % args['indir']):
		if r['gene_id'] != '': genes.add(r['gene_id'])
	return genes

def init_pi(args, samples):
	""" Initialize dictionary to store nucleotide diversity statistics """
	pi = {}
	count_samples = len([_ for _ in samples if _.pass_qc])
	if args['sample_type'] == 'per-sample':
		for s in samples:
			if s.pass_qc:
				if args['genomic_type'] == 'genome-wide':
					pi[s.id] = Diversity()
				else:
					pi[s.id] = {}
					for gene in list_genes(args):
						pi[s.id][gene] = Diversity()
	else:
		if args['genomic_type'] == 'genome-wide':
			pi = Diversity()
			pi.samples = count_samples
		else:
			for gene in list_genes(args):
				pi[gene] = Diversity()
				pi[gene].samples = count_samples
	return pi

def maf(x):
	return x if x <= 0.5 else 1-x

def compute_pi(args, samples):
	pi = init_pi(args, samples)
	index = 0
	for site in diversity.parse_sites(args['indir'], samples): # loop over sites
		
		# stop early
		if index >= args['max_sites']: break
			
		# prune low quality samples for site
		site.prune_samples(args['site_depth'], args['site_ratio'])

		# compute site summary stats
		site.summary_stats(args['weight_by_depth'])
		
		#  filter by prevalence, frequency, site type, # of alt alleles
		filter, code = site.filter_site(args['site_prev'], args['site_maf'], args['site_type'], args['site_freq'])
		if filter: continue
		#  filter 1 - args['rand_sites'] proportion of sites
		elif args['rand_sites'] and random.uniform(0, 1) > args['rand_sites']: continue
		#  keep site
		else: index += 1
			
		# downsample reads
		if args['rand_reads'] and site.maf > 0:
			site.resample_reads(args['rand_reads'], args['replace_reads'])
			site.mean_freq = site.compute_mean_freq(args['weight_by_depth'])
			site.maf = min(site.mean_freq, 1-site.mean_freq)

		# compute pi for pooled-samples
		if args['sample_type'] == 'pooled-samples':
			site_pi = 2 * site.maf * (1-site.maf)
			if args['genomic_type'] == 'genome-wide':
				pi.pi += site_pi
				pi.sites += 1
			else:
				pi[site.info['gene_id']].pi += site_pi
				pi[site.info['gene_id']].sites += 1

		# compute pi per-sample
		else:
			for sample, ref_freq, depth in zip(site.samples, site.ref_freq, site.depth):
				site_maf = maf(ref_freq)
				if args['genomic_type'] == 'genome-wide':
					pi[sample.id].pi += 2 * site_maf * (1-site_maf)
					pi[sample.id].sites += 1
					pi[sample.id].depth += depth
				else:
					pi[sample.id][site.info['gene_id']].pi += 2 * site_maf * (1-site_maf)
					pi[sample.id][site.info['gene_id']].sites += 1
					pi[sample.id][site.info['gene_id']].depth += depth
	return pi

def write_pi(args, samples, pi):
	""" Write nucleotide diversity results to specified output file """
	outfile = open(args['out'], 'w')
	if args['sample_type'] == 'pooled-samples':
		if args['genomic_type'] == 'genome-wide':
			h = ['count_samples', 'count_sites', 'diversity']
			outfile.write('\t'.join([str(_) for _ in h])+'\n')
			r = [pi.samples, pi.sites, pi.pi]
			outfile.write('\t'.join([str(_) for _ in r])+'\n')
		else:
			h = ['gene_id', 'count_samples', 'count_sites', 'diversity']
			outfile.write('\t'.join([str(_) for _ in h])+'\n')
			for gene in pi:
				r = [gene, pi[gene].samples, pi[gene].sites, pi[gene].pi]
				outfile.write('\t'.join([str(_) for _ in r])+'\n')
	elif args['genomic_type'] == 'genome-wide':
		h = ['sample_id', 'depth', 'count_sites', 'diversity']
		outfile.write('\t'.join([str(_) for _ in h])+'\n')
		for s in samples:
			if not s.pass_qc: continue
			r = [s.id, pi[s.id].depth, pi[s.id].sites,  pi[s.id].pi]
			outfile.write('\t'.join([str(_) for _ in r])+'\n')
	else:
		h = ['sample_id', 'gene_id', 'depth', 'count_sites', 'diversity']
		outfile.write('\t'.join([str(_) for _ in h])+'\n')
		for s in samples:
			if not s.pass_qc: continue
			for gene in pi[s.id]:
				r = [s.id, gene, pi[s.id][gene].depth, pi[s.id][gene].sites,  pi[s.id][gene].pi]
				outfile.write('\t'.join([str(_) for _ in r])+'\n')
	outfile.close()

def resample_samples(samples):
	""" Random select high quality samples, set pass_qc=False for others """
	hq_indexes = [index for index, sample in enumerate(samples) if sample.pass_qc]
	if args['rand_samples'] > len(hq_indexes):
		sys.exit("\nError: --rand_samples cannot exceed the number of high-quality samples")
	else:
		hq_indexes = np.random.choice(hq_indexes, args['rand_samples'], replace=False)
		for index, sample in enumerate(samples):
			if index not in hq_indexes: sample.pass_qc = False

def fetch_samples(args):
	""" List and select samples from input """
	from midas.analyze.diversity import Sample
	samples = []
	for info in parse_file('%s/snps_summary.txt' % args['indir']):
		# init sample
		sample = Sample(info)
		sample.pass_qc = True
		# keep/exclude sample
		if sample.filter(args['sample_depth'], args['fract_cov']):
			sample.pass_qc = False
		if args['keep_samples'] and sample.id not in args['keep_samples']:
			sample.pass_qc = False
		if args['exclude_samples'] and sample.id in args['exclude_samples']:
			sample.pass_qc = False
		if sum([1 for s in samples if s.pass_qc]) >= args['max_samples']:
			sample.pass_qc = False
		# store sample
		samples.append(sample)
	# select random sample
	if args['rand_samples']: resample_samples(samples)
	return samples

if __name__ == '__main__':
	args = parse_arguments()
	check_args(args)
	print_copyright()
	print_args(args)
	
	print("\nSelecting subset of samples...")
	samples = fetch_samples(args)

	print("Estimating diversity metrics...")
	pi = compute_pi(args, samples)

	print("Writing results to output file...")
	write_pi(args, samples, pi)


