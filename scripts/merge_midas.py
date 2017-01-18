#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, sys, platform, argparse, numpy as np
from midas import utility

def get_program():
	""" Get program specified by user (species, genes, or snps) """
	if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
		print('Description: merge MIDAS results across metagenomic samples')
		print('')
		print('Usage: merge_midas.py <command> [options]')
		print('')
		print('Commands:')
		print('\tspecies\t build species abundance matrix')
		print('\tgenes\t build pangenome matrix for each species')
		print('\tsnps\t perform multi-sample SNP calling and build SNP matrix for each species')
		print('')
		print('Note: use merge_midas.py <command> -h to view usage for a specific command')
		quit()
	elif sys.argv[1] not in ['species', 'genes', 'snps']:
		sys.exit("\nError: Unrecognized command: '%s'\n" % sys.argv[1])
		quit()
	else:
		return sys.argv[1]

def get_arguments(program):
	""" Get arguments for specified program """
	if program == 'species':
		args = species_arguments()
	elif program == 'genes':
		args = genes_arguments()
	elif program == 'snps':
		args = snps_arguments()
	else:
		sys.exit("\nError: Unrecognized program: '%s'\n" % program)
	return args

def species_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: Merge species abundance files across samples

Usage: merge_midas.py species <outdir> [options]
""",
		epilog="""Examples:
1) provide list of paths to sample directories:
merge_midas.py species /path/to/outdir -i /path/to/samples/sample_1,/path/to/samples/sample_2 -t list

2) provide directory containing all samples:
merge_midas.py species /path/to/outdir -i /path/to/samples -t dir

3) provide file containing paths to sample directories:
merge_midas.py species /path/to/outdir -i /path/to/samples/sample_paths.txt -t file

4) run a quick test:
merge_midas.py species /path/to/outdir -i /path/to/samples -t dir --max_samples 2
	""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str, help='Directory for output files')
	parser.add_argument('-i', type=str, dest='input', required=True,
		help="""Input to sample directories output by run_midas.py; see '-t' for details""")
	parser.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True, metavar="INPUT_TYPE",
		help="""Specify one of the following:
  list: -i is a comma-separated list (ex: /samples/sample_1,/samples/sample_2)
  dir: -i is a directory containing all samples (ex: /samples)
  file: -i is a file of paths to samples (ex: /sample_paths.txt)""")
	parser.add_argument('-d', type=str, dest='db', default=os.environ['MIDAS_DB'] if 'MIDAS_DB' in os.environ else None,
		help="""Path to reference database
By default the MIDAS_DB environmental variable is used""")
	parser.add_argument('--min_cov', metavar='FLOAT', type=float, default=1.0,
		help="""Minimum marker-gene-coverage for estimating species prevalence (1.0)""")
	parser.add_argument('--max_samples', type=int, metavar='INT',
		help="""Maximum number of samples to process.
Useful for testing (use all)""")
	args = vars(parser.parse_args())
	return args

def genes_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: merge results from pan-genome profiling across samples

Usage: merge_midas.py genes <outdir> [options]
""",
		epilog="""Examples:
1) Merge results for all species. Provide list of paths to sample directories:
merge_midas.py genes /path/to/outdir -i sample_1,sample_2 -t list

2) Merge results for one species (id=Bacteroides_vulgatus_57955):
merge_midas.py genes /path/to/outdir --species_id Bacteroides_vulgatus_57955 -i sample_1,sample_2 -t list

3) Exclude low-coverage samples in output matrix:
merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --sample_depth 5.0

4) Use lenient threshold for determining gene presence-absence:
merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --min_copy 0.1

5) Run a quick test:
merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --max_species 1 --max_samples 10
	""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str,
		help="Directory for output files.\nA subdirectory will be created for each species_id")
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""Input to sample directories output by run_midas.py; see '-t' for details""")
	io.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True, metavar="INPUT_TYPE",
		help="""Specify one of the following:
  list: -i is a comma-separated list (ex: /samples/sample_1,/samples/sample_2)
  dir: -i is a directory containing all samples (ex: /samples)
  file: -i is a file of paths to samples (ex: /sample_paths.txt)""")
	io.add_argument('-d', type=str, dest='db', default=os.environ['MIDAS_DB'] if 'MIDAS_DB' in os.environ else None,
		help="""Path to reference database.
By default, the MIDAS_DB environmental variable is used""")
	species = parser.add_argument_group('Species filters (select subset of species from INPUT)')
	species.add_argument('--min_samples', type=int, default=1, metavar='INT',
		help="""All species with >= MIN_SAMPLES (1)""")
	species.add_argument('--species_id', dest='species_id', type=str, metavar='CHAR',
		help="""Comma-separated list of species ids""")
	species.add_argument('--max_species', type=int, metavar='INT',
		help="""Maximum number of species to merge. Useful for testing (use all)""")
	sample = parser.add_argument_group('Sample filters (select subset of samples from INPUT)')
	sample.add_argument('--sample_depth', type=float, default=1.0, metavar='FLOAT',
		help="""Minimum read-depth across all genes with non-zero coverage (1.0)""")
	sample.add_argument('--max_samples', type=int, metavar='INT',
		help="""Maximum number of samples to process. Useful for testing (use all)""")
	gene = parser.add_argument_group('Quantification')
	gene.add_argument('--cluster_pid', type=str, dest='cluster_pid', default='95', choices=['75', '80', '85', '90', '95', '99'],
		help="""In the database, pan-genomes are defined at 6 different %% identity clustering cutoffs.
CLUSTER_PID allows you to quantify gene content for any of these sets of gene clusters.
By default, gene content is reported for genes clustered at 95%% identity
""")
	gene.add_argument('--min_copy', type=float, default=0.35, metavar='FLOAT',
		help="""Genes >= MIN_COPY are classified as present
Genes < MIN_COPY are classified as absent (0.35)""")
	args = vars(parser.parse_args())
	return args

def snps_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: perform multi-sample core-genome SNP calling 

The pipeline can be broken down into the following steps:
  1) take MIDAS output files from multiple samples
  2) identify species to process (based on user criterea, e.g. min # of samples)
  3) scan across the representative genome of each species
  4) pool nucleotide variants from all metagenomic samples & call the major and minor allele
  5) determine if genomic site is a SNP (e.g. minor allele frequency >1%)
  6) determine if genomic site is in the core-genome (e.g. non-zero depth in >95% of samples)
  7) annotate genomic site by gene_id and coding changes
  8) write core-genome SNPs to matrix files

Usage: merge_midas.py snps <outdir> [options]
""",
		epilog="""Examples:
1) Call SNPs for all species. Provide list of paths to sample directories:
merge_midas.py snps /path/to/outdir -i sample_1,sample_2 -t list

2) Call SNPs for one species (id=Bacteroides_vulgatus_57955):
merge_midas.py snps /path/to/outdir --species_id Bacteroides_vulgatus_57955 -i sample_1,sample_2 -t list

3) Merge results for all sites in the core genome, including those that aren't SNPs
(this is useful for comparing core-genome-wide diversity patterns between species):
merge_midas.py snps /path/to/outdir -i /path/to/samples -t dir --core-sites

4) Run a quick test:
merge_midas.py snps /path/to/outdir -i /path/to/samples -t dir --max_species 1 --max_samples 10 --max_sites 1000
	""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str,
		help="Directory for output files. \nA subdirectory will be created for each species_id")
	parser.add_argument('--threads', type=int, default=1, metavar='INT',
		help="Number of CPUs to use (1)")
	parser.add_argument('--sites_per_iter', type=int, default=1000, metavar='INT',
		help=argparse.SUPPRESS)
	parser.add_argument('--max_gb', type=float, metavar='FLOAT',
		help=argparse.SUPPRESS)
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""Input to sample directories output by run_midas.py; see '-t' for details""")
	io.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True, metavar="INPUT_TYPE",
		help="""Specify one of the following:
  list: -i is a comma-separated list (ex: /samples/sample_1,/samples/sample_2)
  dir: -i is a directory containing all samples (ex: /samples)
  file: -i is a file of paths to samples (ex: /sample_paths.txt)""")
	io.add_argument('-d', type=str, dest='db', default=os.environ['MIDAS_DB'] if 'MIDAS_DB' in os.environ else None,
		help="""Path to reference database
By default, the MIDAS_DB environmental variable is used""")
	presets = parser.add_argument_group("Presets (option groups for easily...)")

	snps = parser.add_argument_group("Presets")
	snps.add_argument('--core_snps', action='store_true',
		help="""Same as: --snp_type bi --site_depth 1 --site_ratio 2.0 --site_prev 0.95 (default)""")
	snps.add_argument('--core_sites', action='store_true',
		help="""Same as: --snp_type any --site_depth 1 --site_ratio 2.0 --site_prev 0.95""")
	snps.add_argument('--all_snps', action='store_true',
		help="""Same as: --snp_type bi --site_prev 0.0""")
	snps.add_argument('--all_sites', action='store_true',
		help="""Same as: --snp_type any --site_prev 0.0""")
			
	species = parser.add_argument_group("Species filters (select subset of species from INPUT)")
	species.add_argument('--min_samples', type=int, default=1, metavar='INT',
		help="""All species with >= MIN_SAMPLES (1)""")
	species.add_argument('--species_id', dest='species_id', type=str, metavar='CHAR',
		help="""Comma-separated list of species ids""")
	species.add_argument('--max_species', type=int, metavar='INT',
		help="""Maximum number of species to call SNPs for (all with >= 1 sample)""")
	
	sample = parser.add_argument_group("Sample filters (select subset of samples from INPUT)")
	sample.add_argument('--sample_depth', dest='sample_depth', type=float, default=5.0, metavar='FLOAT',
		help="""Minimum average read depth per sample (5.0)""")
	sample.add_argument('--fract_cov', dest='fract_cov', type=float, default=0.4, metavar='FLOAT',
		help="""Fraction of reference sites covered by at least 1 read (0.4)""")
	sample.add_argument('--max_samples', type=int, metavar='INT',
		help="""Maximum number of samples to process. useful for quick tests (use all)""")
	sample.add_argument('--all_samples', default=False, action='store_true',
		help="""Include all samples in output""")
		
	snps = parser.add_argument_group("Site filters (select subset of genomic sites from INPUT)")
	snps.add_argument('--snp_type', choices=['any', 'mono', 'bi', 'tri', 'quad'], nargs='+', default=['bi'], metavar="",
		help="""Specify one or more of the following:
  mono: keep sites with 1 allele > ALLELE_FREQ
  bi: keep sites with 2 alleles > ALLELE_FREQ (default)
  tri: keep sites with 3 alleles > ALLELE_FREQ
  quad: keep sites with 4 alleles > ALLELE_FREQ
  any: keep sites regardless of observed alleles
""")
	snps.add_argument('--allele_freq', type=float, default=0.01, metavar='FLOAT',
		help="""Minimum frequency for calling an allele present (0.01)
Values > 0.0 and < 0.5 are accepted.
Ex: --snp_type=bi --allele_freq=0.01 keeps bi-allelic SNPs with a minimum frequency of 1%%
""")
	snps.add_argument('--site_depth', type=int, default=1, metavar='INT',
		help="""Minimum number of reads mapped to genomic site (1)
Used in combination with --site_prev to determine if site is in core-genome""")
	snps.add_argument('--site_ratio', type=float, default=2.0, metavar='FLOAT',
		help="""Maximum ratio of site depth to genome depth (2.0)
This filter helps to eliminate genomic sites with abnormally high read depth""")
	snps.add_argument('--site_prev', type=float, default=0.95, metavar='FLOAT',
		help="""Minimum fraction of sample where genomic site is >= SITE_DEPTH and <= SITE_RATIO (0.95)
A high value selects for sites in the core-genome (i.e. present in nearly all strains).
A low value includes sites in variable regions and/or with abnormally high read depth""")
	snps.add_argument('--max_sites', type=int, default=float('Inf'), metavar='INT',
		help="""Maximum number of sites to include in output (use all). Useful for quick tests """)

	args = vars(parser.parse_args())
	args = add_snp_presets(args)
	return args

def add_snp_presets(args):
	""" Set argument values based on selected presets """
	if args['all_samples']:
		args['sample_depth'] = 0.0
		args['fract_cov'] = 0.0
	if args['all_sites']:
		args['site_prev'] = 0.0
		args['snp_type'] = ['any']
	if args['all_snps']:
		args['site_prev'] = 0.0
		args['snp_type'] = ['bi']
	if args['core_sites']:
		args['site_depth'] = 1
		args['site_ratio'] = 2.0
		args['site_prev'] = 0.95
		args['snp_type'] = ['any']
	if args['core_snps']:
		args['site_depth'] = 1
		args['site_ratio'] = 2.0
		args['site_prev'] = 0.95
		args['snp_type'] = ['bi']
	if args['snp_type'] == ['any']:
		args['snp_type'] = ['mono', 'bi', 'tri', 'quad']
	return args

def check_arguments(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program in ['species', 'snps', 'genes']:
		if not os.path.isdir(args['outdir']): os.mkdir(args['outdir'])
		check_input(args)
		utility.check_database(args)
	else:
		sys.exit("\nError: Unrecognized program: '%s'\n" % program)
	if platform.system() not in ['Linux', 'Darwin']:
		sys.exit("\nError: Operating system '%s' not supported\n" % system())

def check_input(args):
	args['indirs'] = []
	error = "\nError: specified input %s does not exist: %s\n"
	if args['intype'] == 'dir':
		if not os.path.isdir(args['input']):
			sys.exit(error % (args['intype'], os.path.abspath(args['input'])))
		else:
			for dir in os.listdir(args['input']):
				args['indirs'].append(os.path.join(args['input'], dir))
	elif args['intype'] == 'file':
		if not os.path.isfile(args['input']):
			sys.exit(error % (args['intype'], os.path.abspath(args['input'])))
		else:
			for line in open(args['input']):
				dir = line.rstrip().rstrip('/')
				if not os.path.isdir(dir): sys.exit(error % ('dir', dir))
				else: args['indirs'].append(dir)
	elif args['intype'] == 'list':
		for dir in args['input'].split(','):
			if not os.path.isdir(dir): sys.exit(error % ('dir', dir))
			else: args['indirs'].append(dir)

def print_arguments(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		print_species_arguments(args)
	elif program == 'genes':
		print_genes_arguments(args)
	elif program == 'snps':
		print_snps_arguments(args)
	else:
		sys.exit("\nError: Unrecognized program: '%s'\n" % program)

def print_species_arguments(args):
	print ("===========Parameters===========")
	print ("Command: %s" % ' '.join(sys.argv))
	print ("Script: merge_midas.py species")
	print ("Database: %s" % args['db'])
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Minimum coverage for estimating prevalence: %s" % args['min_cov'])
	if args['max_samples']: print ("Keep <= %s samples" % args['max_samples'])
	print ("===============================")
	print ("")

def print_genes_arguments(args):
	print ("===========Parameters===========")
	print ("Command: %s" % ' '.join(sys.argv))
	print ("Script: merge_midas.py genes")
	print ("Database: %s" % args['db'])
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Species selection criteria:")
	if args['species_id']: print ("  keep species ids: %s" % args['species_id'].split(','))
	else: print ("  keep species with >= %s samples" % args['min_samples'])
	if args['max_species']: print ("  keep <= %s species" % args['max_species'])
	print ("Sample selection criteria:")
	print ("  keep samples with >=%s mean coverage across genes with non-zero coverage" % args['sample_depth'])
	if args['max_samples']: print ("  keep <= %s samples" % args['max_samples'])
	print ("Gene quantification criterea:")
	print ("  quantify genes clustered at %s%% identity" % args['cluster_pid'])
	print ("  present (1): genes with copy number >= %s" % args['min_copy'])
	print ("  absent (0): genes with copy number < %s" % args['min_copy'])
	print ("===============================")
	print ("")

def print_snps_arguments(args):
	print ("===========Parameters===========")
	print ("Command: %s" % ' '.join(sys.argv))
	print ("Script: merge_midas.py snps")
	print ("Database: %s" % args['db'])
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Species selection criteria:")
	if args['species_id']: print ("  keep species ids: %s" % args['species_id'].split(','))
	else: print ("  keep species with >= %s samples" % args['min_samples'])
	if args['max_species']: print ("  keep <= %s species" % args['max_species'])
	print ("Sample selection criteria:")
	print ("  keep samples with >= %s mean coverage across sites with non-zero coverage" % args['sample_depth'])
	print ("  keep samples where >= %s percent of sites have non-zero coverage" % (100*args['fract_cov']))
	if args['max_samples']: print ("  keep <= %s samples" % args['max_samples'])
	print ("Site selection criteria:")
	print ("  minimum allele frequency for SNP calling: %s%%" % (100*args['allele_freq']))
	print ("  keep %s-alleic sites" % (','.join(args['snp_type'])))
	print ("  keep sites with depth >= %s in >= %s%% of samples" % (args['site_depth'], 100*args['site_prev']) )
	print ("  keep sites with depth <= %sx the mean-genome-wide-depth in >= %s%% of samples" % (args['site_ratio'], 100*args['site_prev']) )
	if args['max_sites'] != float('Inf'): print ("  keep <= %s sites" % (args['max_sites']))
	print ("Number of CPUs to use: %s" % args['threads'])
	print ("===============================")
	print ("")

def run_program(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		from midas.merge import species
		species.run_pipeline(args)
	elif program == 'genes':
		from midas.merge import genes
		genes.run_pipeline(args)
	elif program == 'snps':
		from midas.merge import snps
		snps.run_pipeline(args)
	else:
		sys.exit("\nError: Unrecognized program: '%s'\n" % program)

if __name__ == '__main__':
	program = get_program()
	args = get_arguments(program)
	check_arguments(program, args)
	utility.print_copyright()
	print_arguments(program, args)
	run_program(program, args)



