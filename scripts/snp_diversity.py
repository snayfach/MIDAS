#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, numpy as np, random, csv
from midas.utility import print_copyright
from midas.analyze import parse_snps

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

Usage: snp_diversity.py indir [options]
""",
		epilog="""
Examples:
1) Quantify within-sample heterogenity genome-wide
snp_diversity.py /path/to/snps --genomic_type genome-wide --sample_type per-sample --out /path/to/output

2) Quantify between-sample heterogenity genome-wide
snp_diversity.py /path/to/snps --genomic_type genome-wide --sample_type pooled-sample --out /path/to/output

3) Quantify between-sample heterogenity per-gene
snp_diversity.py /path/to/snps --genomic_type per-gene --sample_type pooled-samples --out /path/to/output

4) Use downsampling to control the read-depth at each genomic site
snp_diversity.py /path/to/snps --genomic_type genome-wide --sample_type per-sample --out /path/to/output

5) Only quantify diversity at non-synonymous sites
snp_diversity.py /path/to/snps --genomic_type genome-wide --sample_type pooled-samples --site_type 4D --locus_type CDS --out /path/to/output

6) Quantify SNPs using a different definition of a polymorphism
snp_diversity.py /path/to/snps --genomic_type genome-wide --sample_type per-sample --snp_maf 0.05 --out /path/to/output

7) Run a quick test
snp_diversity.py /path/to/snps --max_sites 10000  --out /path/to/output

""")
	parser.add_argument('indir', metavar='PATH', type=str,
		help="""path to output from `merge_midas.py snps` for one species
directory should be named according to a species_id and contains files 'snps_*.txt')""")
	parser.add_argument('--out', metavar='PATH', type=str, default='/dev/stdout',
		help="""path to output file (/dev/stdout)""")

	diversity = parser.add_argument_group("Diversity options")
	diversity.add_argument('--genomic_type', choices=['genome-wide', 'per-gene'], default='genome-wide',
		help="""compute diversity for individual genes or genome-wide (genome-wide)""")
	diversity.add_argument('--sample_type', choices=['per-sample', 'pooled-samples'], default='per-sample',
		help="""compute diversity for individual samples or for pooled reads across samples (per-sample)""")
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
	diversity.add_argument('--snp_maf', type=float, metavar='FLOAT', default=0.01,
		help="""minor allele frequency cutoff for determining if a site is a SNP (0.01)""")
	diversity.add_argument('--consensus', action='store_true', default=False,
		help="""call consensus alleles prior to calling SNPs""")

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
	snps.add_argument('--site_list', metavar='PATH',
		help="""path to file containing newline-delimited list of genomic sites to include.\nother filters will still apply""")
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
a value of 10 will filter genomic sites with 10x greater coverage than the genomic background""")
	snps.add_argument('--allele_support', type=float, default=0.5, metavar='FLOAT',
		help="minimum fraction of reads supporting consensus allele (0.50)")
	snps.add_argument('--locus_type', choices=['CDS', 'RNA', 'IGR'],
		help="""use genomic sites that intersect: 'CDS': coding genes, 'RNA': rRNA and tRNA genes, 'IGS': intergenic regions""")
	snps.add_argument('--site_type', choices=['1D','2D','3D','4D'],
		help="""if locus_type == 'CDS', use genomic sites with specified degeneracy: 4D indicates synonymous and 1D non-synonymous sites""")
	snps.add_argument('--max_sites', type=int, default=float('Inf'), metavar='INT',
		help="""maximum number of sites to include in output (use all)
useful for quick tests""")
	args = vars(parser.parse_args())
	if args['rand_reads']: args['site_depth'] = max(args['site_depth'], args['rand_reads'])  
	format_sample_lists(args)
	return args

def format_sample_lists(args):
	keep = args['keep_samples'].rstrip(',').split(',') if args['keep_samples'] else None
	exclude = args['exclude_samples'].rstrip(',').split(',') if args['exclude_samples'] else None

def print_args(args):
	lines = []
	lines.append("Command: %s" % ' '.join(sys.argv))
	lines.append("Script: snp_diversity.py")
	lines.append("Input directory: %s" % args['indir'])
	lines.append("Output file: %s" % args['out'])
	lines.append("Diversity options:")
	lines.append("  genomic_type: %s" % args['genomic_type'])
	lines.append("  sample_type: %s" % args['sample_type'])
	lines.append("  weight_by_depth: %s" % args['weight_by_depth'])
	lines.append("  rand_reads: %s" % args['rand_reads'])
	lines.append("  replace_reads: %s" % args['replace_reads'])
	lines.append("  rand_samples: %s" % args['rand_samples'])
	lines.append("  rand_sites: %s" % args['rand_sites'])
	lines.append("  snp_maf: %s" % args['snp_maf'])
	lines.append("  consensus: %s" % args['consensus'])
	lines.append("Sample filters:")
	lines.append("  sample_depth: %s" % args['sample_depth'])
	lines.append("  fract_cov: %s" % args['fract_cov'])
	lines.append("  max_samples: %s" % args['max_samples'])
	lines.append("  keep_samples: %s" % args['keep_samples'])
	lines.append("  exclude_samples: %s" % args['exclude_samples'])
	lines.append("Site filters:")
	lines.append("  site_list: %s" % args['site_list'])
	lines.append("  site_depth: %s" % args['site_depth'])
	lines.append("  site_prev: %s" % args['site_prev'])
	lines.append("  site_maf: %s" % args['site_maf'])
	lines.append("  site_ratio: %s" % args['site_ratio'])
	lines.append("  allele_support: %s" % args['allele_support'])
	lines.append("  locus_type: %s" % args['locus_type'])
	lines.append("  site_type: %s" % args['locus_type'])	
	lines.append("  max_sites: %s" % args['max_sites'])
	sys.stdout.write('\n'.join(lines)+'\n')

def check_args(args):
	if not os.path.isdir(args['indir']):
		sys.exit("\nError: Specified input directory '%s' does not exist\n" % args['indir'])
	if args['site_depth'] < 2:
		sys.exit("\nError: --site_depth must be >=2 to calculate nucleotide variation\n")
	if args['max_sites'] < 1:
		sys.exit("\nError: --max_sites must be >= 1 to calculate nucleotide variation\n")
	if args['max_samples'] < 1:
		sys.exit("\nError: --max_samples must be >= 1 to calculate nucleotide variation\n")
	if args['site_ratio'] < 0:
		sys.exit("\nError: --site_ratio cannot be a negative number\n")
	if args['site_depth'] < 0:
		sys.exit("\nError: --site_depth cannot be a negative number\n")
	if args['sample_depth'] < 0:
		sys.exit("\nError: --sample_depth cannot be a negative number\n")
	if not 0 <= args['site_maf'] <= 1:
		sys.exit("\nError: --site_maf must be between 0 and 1\n")
	if not 0 <= args['site_prev'] <= 1:
		sys.exit("\nError: --site_prev must be between 0 and 1\n")
	if not 0 <= args['fract_cov'] <= 1:
		sys.exit("\nError: --fract_cov must be between 0 and 1\n")
	if args['rand_reads'] > args['site_depth'] and not args['replace_reads']:
		sys.exit("\nError: --rand_reads cannot exceed --site_depth when --replace_reads=False\n")
	if args['rand_sites'] and (args['rand_sites'] < 0 or args['rand_sites'] > 1):
		sys.exit("\nError: --rand_sites must be between 0 and 1\n")
	if args['locus_type'] != 'CDS' and args['genomic_type'] == 'per-gene':
		sys.exit("\nError: --locus_type must be CDS if --genomic_type is per-gene\n")
	if args['locus_type'] != 'CDS' and args['site_type'] is not None:
		sys.exit("\nError: --locus_type must be CDS if --site_type is specified\n")

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
	path = '%s/snps_info.txt' % args['indir']
	for r in csv.DictReader(open(path), delimiter='\t'):
		if r['gene_id'] != '': genes.add(r['gene_id'])
	return genes

def init_pi(args, samples):
	""" Initialize dictionary to store nucleotide diversity statistics """
	pi = {}
	if args['sample_type'] == 'per-sample':
		if args['genomic_type'] == 'genome-wide':
			for s in samples.values(): pi[s.id] = Diversity()
		elif args['genomic_type'] == 'per-gene':
			genes = list_genes(args)
			for s in samples.values():
				pi[s.id] = {}
				for gene in genes: pi[s.id][gene] = Diversity()
	elif args['sample_type'] == 'pooled-samples':
		if args['genomic_type'] == 'genome-wide':
			pi = Diversity()
			pi.samples = len(samples)
		elif args['genomic_type'] == 'per-gene':
			genes = list_genes(args)
			for gene in genes:
				pi[gene] = Diversity()
				pi[gene].samples = len(samples)
	return pi

def compute_maf(freq):
	""" Compute minor allele frequency """
	return min(freq, 1-freq)

def compute_pi(freq):
	""" Compute diversity based on minor allele frequency """
	return 2*freq*(1-freq)

def is_snp(freq, min_maf):
	""" Determine if a genomic site is a SNP or not """
	maf = compute_maf(freq)
	if maf >= min_maf:
		return True
	else:
		return False

def compute_snp_diversity(args, species, samples, progress):

	pi = init_pi(args, samples)
	
	# read list of genomic sites to keep
	if args['site_list']:
		site_list = [_.rstrip() for _ in open(args['site_list'])]
		site_index = 0
	
	index = 0
	for site in parse_snps.fetch_sites(species, samples):
				
		# print progress			
		if progress and not int(site.id) % 10000:
			sys.stdout.write(" %s sites processed\r" % site.id)
			sys.stdout.flush()

		# stop early
		if index >= args['max_sites']: break
		
		# skip sites not in site_list
		if args['site_list']:
			if site_index >= len(site_list):
				break
			elif site.id != site_list[site_index]:
				continue
			else:
				site_index += 1
			
		#  skip random subset of genomic sites
		if args['rand_sites'] and random.uniform(0, 1) > args['rand_sites']:
			continue
			
		# prune low quality samples for site:
		#   site.samples['sample'].keep = [True|False]
		site.flag_samples(args['site_depth'], args['site_ratio'], args['allele_support'])
						
		# call consensus
		if args['consensus']:
			site.call_consensus()
		
		# compute site summary stats
		#   site.prevalence
		#   site.pooled_maf
		site.summary_stats(args['weight_by_depth'])
		
		# filter genomic site
		#   site.keep = [True|False]
		site.filter(args['site_prev'], args['site_maf'], args['locus_type'], args['site_type'])
		if not site.keep:
			continue
		else:
			index += 1

		# downsample reads & recompute pooled frequency
		if args['rand_reads'] and site.pooled_maf > 0.0:
			site.resample_reads(args['rand_reads'], args['replace_reads'])
			site.pooled_maf = site.compute_pooled_maf(args['weight_by_depth'])

		# compute pi for pooled-samples
		if args['sample_type'] == 'pooled-samples':
			if args['genomic_type'] == 'genome-wide':
				pi.pi += compute_pi(site.pooled_maf)
				pi.snps += 1 if is_snp(site.pooled_maf, args['snp_maf']) else 0
				pi.sites += 1
			else:
				pi[site.gene_id].pi += compute_pi(site.pooled_maf)
				pi[site.gene_id].snps += 1 if is_snp(site.pooled_maf, args['snp_maf']) else 0
				pi[site.gene_id].sites += 1
	
		# compute pi per-sample
		else:
			for sample in site.samples.values():
				if sample.keep:
					if args['genomic_type'] == 'genome-wide':			
						pi[sample.id].pi += compute_pi(sample.freq)
						pi[sample.id].snps += 1 if is_snp(sample.freq, args['snp_maf']) else 0
						pi[sample.id].sites += 1
						pi[sample.id].depth += sample.depth
					else:
						pi[sample.id][site.gene_id].pi += compute_pi(sample.freq)
						pi[sample.id][site.gene_id].snps += 1 if is_snp(sample.freq, args['snp_maf']) else 0
						pi[sample.id][site.gene_id].sites += 1
						pi[sample.id][site.gene_id].depth += sample.depth
								
	return pi

def write_pi(args, samples, pi):
	""" Write nucleotide diversity results to specified output file """
	outfile = open(args['out'], 'w')
	if args['sample_type'] == 'pooled-samples':
		if args['genomic_type'] == 'genome-wide':
			h = ['samples', 'sites', 'snps', 'pi', 'snps_kb', 'pi_bp']
			outfile.write('\t'.join([str(_) for _ in h])+'\n')
			snps_kb = 1000*pi.snps/float(pi.sites) if pi.sites > 0 else 'NA'
			pi_bp = pi.pi/float(pi.sites) if pi.sites > 0 else 'NA'
			r = [pi.samples, pi.sites, pi.snps, pi.pi, snps_kb, pi_bp]
			outfile.write('\t'.join([str(_) for _ in r])+'\n')
		else:
			h = ['gene_id', 'samples', 'sites', 'snps', 'pi', 'snps_kb', 'pi_bp']
			outfile.write('\t'.join([str(_) for _ in h])+'\n')
			for gene in pi:
				snps_kb = 1000*pi[gene].snps/float(pi[gene].sites) if pi[gene].sites > 0 else 'NA'
				pi_bp = pi[gene].pi/float(pi[gene].sites) if pi[gene].sites > 0 else 'NA'
				r = [gene, pi[gene].samples, pi[gene].sites, pi[gene].snps, pi[gene].pi, snps_kb, pi_bp]
				outfile.write('\t'.join([str(_) for _ in r])+'\n')
	elif args['genomic_type'] == 'genome-wide':
		h = ['sample_id', 'depth', 'sites', 'snps', 'pi', 'snps_kb', 'pi_bp']
		outfile.write('\t'.join([str(_) for _ in h])+'\n')
		for s in samples.values():
			snps_kb = 1000*pi[s.id].snps/float(pi[s.id].sites) if pi[s.id].sites > 0 else 'NA'
			pi_bp = pi[s.id].pi/float(pi[s.id].sites) if pi[s.id].sites > 0 else 'NA'
			r = [s.id, pi[s.id].depth, pi[s.id].sites,  pi[s.id].snps, pi[s.id].pi, snps_kb, pi_bp]
			outfile.write('\t'.join([str(_) for _ in r])+'\n')
	else:
		h = ['sample_id', 'gene_id', 'depth', 'sites', 'snps', 'pi', 'snps_kb', 'pi_bp']
		outfile.write('\t'.join([str(_) for _ in h])+'\n')
		for s in samples.values():
			for gene in pi[s.id]:
				snps_kb = 1000*pi[s.id][gene].snps/float(pi[s.id][gene].sites) if pi[s.id][gene].sites > 0 else 'NA'
				pi_bp = pi[s.id][gene].pi/float(pi[s.id][gene].sites) if pi[s.id][gene].sites > 0 else 'NA'
				r = [s.id, gene, pi[s.id][gene].depth, pi[s.id][gene].sites,  pi[s.id][gene].snps, pi[s.id][gene].pi, snps_kb, pi_bp]
				outfile.write('\t'.join([str(_) for _ in r])+'\n')
	outfile.close()

if __name__ == '__main__':
	args = parse_arguments()
	check_args(args)
	print_copyright()
	print_args(args)

	print("\nSelecting subset of samples...")
	species = parse_snps.Species(args['indir'])
	samples = parse_snps.fetch_samples(species, args['sample_depth'], args['fract_cov'], args['max_samples'],
						    		   args['keep_samples'], args['exclude_samples'], args['rand_samples'])
	print(" %s samples selected" % len(samples))

	print("Estimating diversity metrics...\n")
	pi = compute_snp_diversity(args, species, samples, progress=False)

	write_pi(args, samples, pi)


