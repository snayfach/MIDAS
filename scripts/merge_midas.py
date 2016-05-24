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
		sys.error("Unrecognized program: '%s'" % program)
	args = utility.add_ref_db(args)
	return args

def species_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: merge_midas.py species [options]

Description: Merge species abundance files across samples
Input: list of sample directories
Output: relative abundance matrix, genome-coverage matrix, read-count matrix, species prevalence
""",
		epilog="""Examples:
1) provide list of paths to sample directories:
merge_midas.py species  -i /path/to/samples/sample_1,/path/to/samples/sample_2 -t list -o outbase

2) provide directory containing all samples:
merge_midas.py species  -i /path/to/samples -t dir -o outbase

3) provide file containing paths to sample directoriess:
merge_midas.py species  -i /path/to/samples/sample_paths.txt -t file -o outbase
""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('-i', type=str, dest='input', required=True,
		help="""input to sample directories output by run_midas.py species
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
	parser.add_argument('-o', dest='outbase', type=str, required=True,
		help="basename for output files")
	parser.add_argument('-m', dest='min_cov', type=float, required=False, default=1.0,
		help="""minimum genome-coverage for estimating species prevalence (1.0)""")
	args = vars(parser.parse_args())
	return args

def genes_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: merge_midas.py genes [options]

Description: merge results from pan-genome profiling across samples
Input: list of sample directories
Output: pan-genome copy-number matrix, presence/absence matrix, and read-depth matrix
        matrixes also created for KEGG, FIGfams, Gene Ontology, and Enzyme Comission (E.C.)
""",
		epilog="""Examples:
1) Merge results for species 57955. Provide list of paths to sample directories:
merge_midas.py genes -s 57955 -o outdir/57955 -i sample_1,sample_2 -t list

2) Build matrix for pan-genome genes at lower percent id threshold:
merge_midas.py genes -s 57955 -o outdir/57955 -i /path/to/samples -t dir --cluster_pid 85

3) Exclude low-coverage samples in output matrix:
merge_midas.py genes -s 57955 -o outdir/57955 -i /path/to/samples -t dir --marker_coverage 5.0

4) Use lenient threshold for determining gene presence-absence:
merge_midas.py genes -s 57955 -o outdir/57955 -i /path/to/samples -t dir --min_copy 0.1

""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""input to sample directories output by run_midas.py genes
see '-t' for details""")
	io.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True,
		help="""'list': -i is a comma-separated list of paths to sample directories (ex: /sample1,/sample2)
'dir': -i is a  directory containing all samples (ex: /samples_dir)
'file': -i is a file containing paths to sample directories (ex: sample_paths.txt)
""")
	io.add_argument('-o', type=str, dest='outdir', required=True,
		help="""output directory""")
	species = parser.add_argument_group('Species filters (select subset of species from INPUT)')
	species.add_argument('--min_samples', type=int, default=1, metavar='INT',
		help="""all species with >= MIN_SAMPLES (1)""")
	species.add_argument('--sp_id', dest='species_id', type=str,
		help="""comma-separated list of species ids
a list of prevalent species can be obtained by running 'merge_midas.py species'
a map of species ids to species names can be found in 'ref_db/annotations.txt'""")
	species.add_argument('--max_species', type=int, metavar='INT',
		help="""maximum number of species to merge. useful for testing (use all)""")
	sample = parser.add_argument_group('Sample filters (select subset of samples from INPUT)')
	#sample.add_argument('--marker_coverage', type=float, default=1.0, metavar='FLOAT',
	#	help="""minimum coverage per sample across 15 phylogenetic marker genes (1.0)""")
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
Usage: merge_midas.py snps [options]

Description: merge single-nucleotide variants for an individual species across samples
Input: list of sample directories
Output: core-genome SNPs, SNP annotations, SNP allele frequency matrix, SNP alternate alleles, SNP depth,
        core-genome consensus sequences, and a phylogenetic tree
""",
		epilog="""Examples:
1) Merge results for species 57955. Provide list of paths to sample directories:
merge_midas.py snps -s 57955 -o outdir/57955 -i sample_1,sample_2 -t list

2) Only use samples with >15x average depth and only use sites covered by >=10 reads in at least >=95% of samples:
merge_midas.py snps -s 57955 -o outdir/57955 -i /path/to/samples -t dir --sample_depth 15 --site_depth 10 --site_prev 95

3) Just identify core-genome sites and build matrixes; do not build consensus seqs or phylogenetic trees:
merge_midas.py snps -s 57955 -o outdir/57955 -i /path/to/samples -t dir --snps --freq

4) Run a quick test:
merge_midas.py snps -s 57955 -o outdir/57955 -i /path/to/samples -t dir --max_samples 10 --max_sites 1000
""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('--threads', type=int, default=1, metavar='INT',
		help="number of CPUs to use for merging files (1)\nincreases speed when merging across many samples")
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""input to sample directories output by run_midas.py genes
see '-t' for details""")
	io.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True,
		help="""'list': -i is a comma-separated list of paths to sample directories (ex: /sample1,/sample2)
'dir': -i is a  directory containing all samples (ex: /samples_dir)
'file': -i is a file containing paths to sample directories (ex: sample_paths.txt)
""")
	io.add_argument('-o', type=str, dest='outdir', required=True,
		help="""output directory""")
	species = parser.add_argument_group("Species filters (select subset of species from INPUT)")
	species.add_argument('--min_samples', type=int, default=1, metavar='INT',
		help="""all species with >= MIN_SAMPLES (1)""")
	species.add_argument('--sp_id', dest='species_id', type=str,
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
		sys.error("Unrecognized program: '%s'" % program)
	if platform.system() not in ['Linux', 'Darwin']:
		sys.exit("Operating system '%s' not supported" % system())

def check_species(args):
	base_dir = os.path.dirname(os.path.abspath(args['outbase']))
	if not os.path.isdir(base_dir):
		sys.exit("Output directory of base name (-o) not found:\n%s" % base_dir)

def check_genes(args):
	if not os.path.isdir(args['outdir']): os.mkdir(args['outdir'])
	check_input(args)

def check_snps(args):
	if not os.path.isdir(args['outdir']):
		os.mkdir(args['outdir'])
#	if not any([args['snps'], args['freq'], args['cons'], args['tree']]):
#		args['snps'] = True
#		args['freq'] = True
#		args['cons'] = True
#		args['tree'] = True
	check_input(args)

def check_input(args):
	error = "\nError: specified input does not exist: %s"
	if (args['intype'] == 'dir'
			and not os.path.isdir(args['input'])):
		sys.exit(error % os.path.abspath(args['input']))
	elif (args['intype'] == 'file'
			and not os.path.isfile(args['input'])):
		sys.exit(error % os.path.abspath(args['input']))
	elif args['intype'] == 'list':
		for file in input.split(','):
			if not os.path.isfile(file):
				sys.exit(error % file)

def print_arguments(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		print_species_arguments(args)
	elif program == 'genes':
		print_genes_arguments(args)
	elif program == 'snps':
		print_snps_arguments(args)
	else:
		sys.error("Unrecognized program: '%s'" % program)

def print_species_arguments(args):
	print ("===========Parameters===========")
	print ("Script: merge_midas.py species")
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output basename: %s" % args['outbase'])
	print ("Minimum coverage for estimating prevalence: %s" % args['min_cov'])
	print ("")

def print_genes_arguments(args):
	print ("===========Parameters===========")
	print ("Script: merge_midas.py genes")
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Species selection criteria:")
	if args['species_id']:
		print ("  species_ids: %s" % args['species_id'].split(','))
	if args['min_samples']:
		print ("  >= %s high-coverage samples per species" % args['min_samples'])
	if args['max_species']:
		print ("  analyze up to %s species" % args['max_species'])
	print ("Sample selection criteria:")
	#if args['marker_coverage']:
	#	print ("  >=%s average coverage across 15 universal-single-copy genes" % args['marker_coverage'])
	if args['sample_depth']:
		print ("  >=%s average read depth across detected genes" % args['sample_depth'])
	if args['max_samples']:
		print ("  analyze up to %s samples" % args['max_samples'])
	print ("Gene quantification criterea:")
	print ("  present (1): genes with copy number >=%s" % args['min_copy'])
	print ("  absent (0): genes with copy number <%s" % args['min_copy'])
	print ("  cluster genes at %s percent identity" % args['cluster_pid'])
	print ("")

def print_snps_arguments(args):
	print ("===========Parameters===========")
	print ("Script: merge_midas.py snps")
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Species identifier: %s" % args['species_id'])
	print ("Number of CPUs to use: %s" % args['threads'])
	print ("Sample selection criteria:")
	if args['sample_depth']:
		print ("  keep samples with >=%s average coverage across reference genome" % args['sample_depth'])
	if args['fract_cov']:
		print ("  keep samples where >=%s percent of reference genome has non-zero coverage" % (100*args['fract_cov']))
	if args['max_samples']:
		print ("  analyze up to %s samples" % args['max_samples'])
	print ("Site selection criteria:")
	print ("  site must be covered by at least %s reads across %s percent of samples" % (args['site_depth'], 100*args['site_prev']))
	if args['max_sites'] != float('Inf'):
		print ("  analyze up to %s sites" % (args['max_sites']))
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
		sys.error("Unrecognized program: '%s'" % program)

if __name__ == '__main__':
	program = get_program()
	args = get_arguments(program)
	check_arguments(program, args)
	utility.print_copyright()
	print_arguments(program, args)
	run_program(program, args)



