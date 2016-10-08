#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, sys, platform, argparse, numpy as np
from midas import utility

def get_program():
	""" Get program specified by user (species, genes, or snps) """
	if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
		print('')
		print('Usage: merge_midas.py <command> [options]')
		print('')
		print('Note: use merge_midas.py <command> -h to view usage for a specific command')
		print('')
		print('Commands:')
		print('\tspecies\t merge abundances of bacterial species across samples')
		print('\tgenes\t merge pan-genome gene copy numbers of species across samples')
		print('\tsnps\t merge single nucleotide variants of species across samples')
		quit()
	elif sys.argv[1] not in ['species', 'genes', 'snps']:
		sys.exit("Unrecognized command: '%s'" % sys.argv[1])
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
		sys.exit("Unrecognized program: '%s'" % program)
	if 'db' not in args or args['db'] is None:
		utility.add_ref_db(args)
	return args

def species_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: Merge species abundance files across samples

Usage: merge_midas.py species outdir [options]
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

Output files:
1) relative_abundance.txt: relative abundance matrix (columns are samples, rows are species)
2) count_reads.txt: read count matrix (columns are samples, rows are species)
3) coverage.txt: genome coverage matrix (columns are samples, rows are species)
4) species_prevalence.txt: summary statistics for each species across samples

Output formats:
species_prevalence.txt
1) species_id: species identifier
2) species_name: unique species name
3) mean_coverage: average read-depth across samples
4) median_coverage: median read-depth across samples
5) mean_abundance: average relative abundance across samples
6) median_abundance: median relative abundance across samples
7) prevalence: number of samples with >= `MIN_COV`
""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str, help='directory for output files')
	parser.add_argument('-i', type=str, dest='input', required=True,
		help="""input to sample directories output by run_midas.py
can be a list of directories, a directory containing all samples, or a file with paths
see '-t' for details""")
	parser.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True,
		help="""'list': -i incdicates a comma-separated list of paths to sample directories
        example: /path/to/samples/sample_1,/path/to/samples/sample_2
'dir': -i incdicates a  directory containing all samples
       example: /path/to/samples
'file': -i incdicates a file containing paths to sample directories
	   example: /path/to/sample_paths.txt
""")
	parser.add_argument('--ref_db', dest='db', type=str, metavar='PATH',
		help="""path to alternative reference database directory""")
	parser.add_argument('--min_cov', metavar='FLOAT', type=float, default=1.0,
		help="""minimum marker-gene-coverage for estimating species prevalence (1.0)""")
	parser.add_argument('--max_samples', type=int, metavar='INT',
		help="""maximum number of samples to process. useful for testing (use all)""")
	args = vars(parser.parse_args())
	return args

def genes_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: merge results from pan-genome profiling across samples

Usage: merge_midas.py genes outdir [options]
""",
		epilog="""Examples:
1) Merge results for all species. Provide list of paths to sample directories:
merge_midas.py genes /path/to/outdir -i sample_1,sample_2 -t list

2) Merge results for one species (id=57955):
merge_midas.py genes /path/to/outdir --species_id 57955 -i sample_1,sample_2 -t list

3) Build matrix for pan-genome genes at lower percent id threshold:
merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --cluster_pid 85

4) Exclude low-coverage samples in output matrix:
merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --sample_depth 5.0

5) Use lenient threshold for determining gene presence-absence:
merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --min_copy 0.1

6) Run a quick test:
merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --max_species 1 --max_samples 10

Output files:
1) genes_copynum.txt: gene copy-number matrix (columns are samples, rows are gene families)
2) genes_presabs.txt: gene presence-absence (0/1) matrix (columns are samples, rows are gene families). genes with copy-number >= MIN_COPY are called as present and genes below MIN_COPY are called as absent
3) genes_depth.txt: gene coverage (i.e. read depth) matrix (columns are samples, rows are gene families)
4) genes_info.txt: detailed information of genes
5) genes_summary.txt: alignment summary statistics of genes across samples

Output formats:
genes_info.txt
1) gene_id: identifier of 99% identity gene family
2) family_id: mapping to gene family clustered at CLUSTER_PID
3) function_id: identifier of function
4) function_db: database (kegg, figfam, go, ec) corresponding to function_id

genes_summary.txt
1) sample_id: sample identifier
2) pangenome_size: total number of gene families (99% identity clustering cutoff) in reference pangenome
3) covered_genes: number of pangenome gene families with non-zero coverage
4) fraction_covered: fraction of pangenome gene families with non-zero coverage
5) mean_coverage: mean read-depth across gene families with non-zero coverage
6) marker_coverage: median read-depth across 15 universal single copy genes

""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str,
		help="directory for output files. a subdirectory will be created for each species_id")
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""input to sample directories output by run_midas.py
see '-t' for details""")
	io.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True,
		help="""'list': -i is a comma-separated list of paths to sample directories (ex: /sample1,/sample2)
'dir': -i is a  directory containing all samples (ex: /samples_dir)
'file': -i is a file containing paths to sample directories (ex: sample_paths.txt)
""")
	io.add_argument('--ref_db', dest='db', type=str, metavar='PATH',
		help="""path to alternative reference database directory""")
	species = parser.add_argument_group('Species filters (select subset of species from INPUT)')
	species.add_argument('--min_samples', type=int, default=1, metavar='INT',
		help="""all species with >= MIN_SAMPLES (1)""")
	species.add_argument('--species_id', dest='species_id', type=str, metavar='CHAR',
		help="""comma-separated list of species ids
a list of prevalent species can be obtained by running 'merge_midas.py species'
a map of species ids to species names can be found in 'ref_db/annotations.txt'""")
	species.add_argument('--max_species', type=int, metavar='INT',
		help="""maximum number of species to merge. useful for testing (use all)""")
	sample = parser.add_argument_group('Sample filters (select subset of samples from INPUT)')
	sample.add_argument('--sample_depth', type=float, default=1.0, metavar='FLOAT',
		help="""minimum coverage per sample across all genes with non-zero coverage (1.0)""")
	sample.add_argument('--max_samples', type=int, metavar='INT',
		help="""maximum number of samples to process. useful for testing (use all)""")
	gene = parser.add_argument_group('Presence/Absence')
	gene.add_argument('--min_copy', type=float, default=0.35, metavar='FLOAT',
		help="""genes >= MIN_COPY are classified as present
genes < MIN_COPY are classified as absent (0.35)""")
	gene.add_argument('--cluster_pid', type=str, dest='cluster_pid', default='95', choices=['75', '80', '85', '90', '95', '99'],
		help="""gene family percent identity
small values: fewer, larger gene families
large values: more, smaller gene families (95)""")
	args = vars(parser.parse_args())
	return args

def snps_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: merge single-nucleotide variant results across samples

Usage: merge_midas.py snps outdir [options]
""",
		epilog="""Examples:
1) Merge results for all species. Provide list of paths to sample directories:
merge_midas.py snps /path/to/outdir -i sample_1,sample_2 -t list

2) Merge results for one species (id=57955):
merge_midas.py snps /path/to/outdir --species_id 57955 -i sample_1,sample_2 -t list

3) Only use samples with >15x average depth and only use sites covered by >=10 reads in at least >=95% of samples:
merge_midas.py snps /path/to/outdir -i /path/to/samples -t dir --sample_depth 15 --site_depth 10 --site_prev 0.95

4) Run a quick test:
merge_midas.py snps /path/to/outdir -i /path/to/samples -t dir --max_species 1 --max_samples 10 --max_sites 1000

Output files:
1) snps_ref_freq.txt: reference allele frequency matrix (sites x samples). each value is the proportion of reads that matched the reference allele for a sample at a genomic site
2) snps_alt_allele.txt: alternate allele matrix (sites x samples). each value is the alternate allele observed for a sample at a genomic site
3) snps_depth.txt: site depth matrix (sites x samples). each value is the total number of reads observed for a sample at a genomic site
4) snps_info.txt: detailed information for each genomic site included in output
5) snps_summary.txt: alignment summary statistics for all samples

File formats:
snps_info.txt
1) site_id: identifier of genomic site
2) mean_freq: average frequency of reference allele across samples
3) mean_depth: average number of mapped reads across samples
4) site_prev: proportion of samples with sufficient depth
5) ref_allele: reference allele
6) allele_props: distribution of allele frequencies across samples
7) site_type: NC=non-coding site, 1D-4D=coding site
8) gene_id: gene identifier
9) amino_acids: amino acid for each possible allele
10) snps: indicates whether an allele is synonymous (SYN) or non-synonymous (NS)

snps_summary.txt
1) sample_id: sample identifier
2) genome_length: length of reference genome used for read-mapping
3) covered_bases: number of genomic positions covered
4) fraction_covered: fraction of genomic positions with non-zero coverage
5) mean_coverage: mean read-depth at covered genomic positions

""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str,
		help="directory for output files. a subdirectory will be created for each species_id")
	parser.add_argument('--threads', type=int, default=1, metavar='INT',
		help="number of CPUs to use for merging files (1)\nincreases speed when merging across many samples")
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""input to sample directories output by run_midas.py
see '-t' for details""")
	io.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True,
		help="""'list': -i is a comma-separated list of paths to sample directories (ex: /sample1,/sample2)
'dir': -i is a  directory containing all samples (ex: /samples_dir)
'file': -i is a file containing paths to sample directories (ex: sample_paths.txt)
""")
	io.add_argument('--ref_db', dest='db', type=str, metavar='PATH',
		help="""path to alternative reference database directory""")
	species = parser.add_argument_group("Species filters (select subset of species from INPUT)")
	species.add_argument('--min_samples', type=int, default=1, metavar='INT',
		help="""all species with >= MIN_SAMPLES (1)""")
	species.add_argument('--species_id', dest='species_id', type=str, metavar='CHAR',
		help="""comma-separated list of species ids
a list of prevalent species can be obtained by running 'merge_midas.py species'
a map of species ids to species names can be found in 'ref_db/annotations.txt'""")
	species.add_argument('--max_species', type=int, metavar='INT',
		help="""maximum number of species to merge (use all)""")
	sample = parser.add_argument_group("Sample filters (select subset of samples from INPUT)")
	sample.add_argument('--sample_depth', dest='sample_depth', type=float, default=5.0, metavar='FLOAT',
		help="""minimum average read depth per sample (5.0)""")
	sample.add_argument('--fract_cov', dest='fract_cov', type=float, default=0.4, metavar='FLOAT',
		help="""fraction of reference sites covered by at least 1 read (0.4)""")
	sample.add_argument('--max_samples', type=int, metavar='INT',
		help="""maximum number of samples to process.
useful for quick tests (use all)""")
	snps = parser.add_argument_group("Site filters (select subset of genomic sites from INPUT)")
	snps.add_argument('--site_depth', type=int, default=3, metavar='INT',
		help="""minimum number of mapped reads per site.
a high value like 20 will result in accurate allele frequencies, but may discard many sites.
a low value like 1 will retain many sites but may not result in accurate allele frequencies (3)""")
	snps.add_argument('--site_prev', type=float, default=0.95, metavar='FLOAT',
		help="""site has at least <site_depth> coverage in at least <site_prev> proportion of samples.
a value of 1.0 will select sites that have sufficent coverage in all samples.
a value of 0.0 will select all sites, including those with low coverage in many samples 
NAs recorded for included sites with less than <site_depth> in a sample (0.95)""")
	snps.add_argument('--site_maf', type=float, default=0.0, metavar='FLOAT',
		help="""minimum minor allele frequency of site across samples.
setting this to zero (default) will keep invariant sites across samples.
setting this above zero (e.g. 0.01, 0.02, 0.05) will only keep common variants""")
	snps.add_argument('--max_sites', type=int, default=float('Inf'), metavar='INT',
		help="""maximum number of sites to include in output.
useful for quick tests (use all)""")
	args = vars(parser.parse_args())
	return args

def check_arguments(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		check_species(args)
	elif program == 'genes':
		check_genes(args)
	elif program == 'snps':
		check_snps(args)
	else:
		sys.exit("Unrecognized program: '%s'" % program)
	if platform.system() not in ['Linux', 'Darwin']:
		sys.exit("Operating system '%s' not supported" % system())

def check_species(args):
	if not os.path.isdir(args['outdir']): os.mkdir(args['outdir'])
	check_input(args)

def check_genes(args):
	if not os.path.isdir(args['outdir']): os.mkdir(args['outdir'])
	check_input(args)

def check_snps(args):
	if not os.path.isdir(args['outdir']):
		os.mkdir(args['outdir'])
	check_input(args)

def check_input(args):
	args['indirs'] = []
	error = "\nError: specified input %s does not exist: %s"
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
		sys.exit("Unrecognized program: '%s'" % program)

def print_species_arguments(args):
	print ("===========Parameters===========")
	print ("Script: merge_midas.py species")
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Minimum coverage for estimating prevalence: %s" % args['min_cov'])
	if args['max_samples']: print ("Keep <= %s samples" % args['max_samples'])
	print ("")

def print_genes_arguments(args):
	print ("===========Parameters===========")
	print ("Script: merge_midas.py genes")
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
	print ("  present (1): genes with copy number >= %s" % args['min_copy'])
	print ("  absent (0): genes with copy number < %s" % args['min_copy'])
	print ("  cluster genes at %s percent identity" % args['cluster_pid'])
	print ("")

def print_snps_arguments(args):
	print ("===========Parameters===========")
	print ("Script: merge_midas.py snps")
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
	print ("  keep sites covered by >= %s reads across >= %s percent of samples" % (args['site_depth'], 100*args['site_prev']))
	if args['max_sites'] != float('Inf'): print ("  keep <= %s sites" % (args['max_sites']))
	print ("Number of CPUs to use: %s" % args['threads'])
	print ("")

def run_program(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		from midas.merge import merge_species
		merge_species.run_pipeline(args)
	elif program == 'genes':
		from midas.merge import merge_genes
		merge_genes.run_pipeline(args)
	elif program == 'snps':
		from midas.merge import merge_snps
		merge_snps.run_pipeline(args)
	else:
		sys.exit("Unrecognized program: '%s'" % program)

if __name__ == '__main__':
	program = get_program()
	args = get_arguments(program)
	check_arguments(program, args)
	utility.print_copyright()
	print_arguments(program, args)
	run_program(program, args)



