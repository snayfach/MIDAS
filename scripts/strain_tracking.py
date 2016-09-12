#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, platform, itertools
from operator import itemgetter
from midas import utility

def get_program():
	""" Get program specified by user (species, genes, or snps) """
	if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
		print('')
		print('Usage: strain_tracking.py <command> [options]')
		print('')
		print('Note: use strain_tracking.py <command> -h to view usage for a specific command')
		print('')
		print('Commands:')
		print('\tid_markers      identify rare SNPs that disriminate individual strains')
		print('\ttrack_markers   track rare SNPs between samples and determine transmission')
		quit()
	elif sys.argv[1] not in ['id_markers', 'track_markers']:
		sys.exit("Unrecognized command: '%s'" % sys.argv[1])
		quit()
	else:
		return sys.argv[1]

def get_arguments(program):
	""" Get arguments for specified program """
	if program == 'id_markers':
		args = id_arguments()
	elif program == 'track_markers':
		args = track_arguments()
	else:
		sys.exit("Unrecognized program: '%s'" % program)
	return args

def id_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: identify rare SNPs that disriminate individual strains of a particular species

Usage: strain_tracking.py id_markers [options]
""",
	epilog="""
Examples:
1) Identify marker alleles for species_id 57955 which are found in exactly one sample
strain_tracking.py id_markers --indir merged_snps/57955 --out 57955.markers --max_groups 1

2) Run a quick test for species_id 57955
strain_tracking.py id_markers --indir merged_snps/57955 --out 57955.markers --max_sites 10000

Output fields:
  site_id: site identifier
  allele: nucleotide (A, T, C, or G)
  count_samples - number of sample-groups with non-zero coverage at site_id
  count_A: number of sample-groups with A at site_id
  count_T: number of sample-groups with T at site_id
  count_C: number of sample-groups with C at site_id
  count_G: number of sample-groups with G at site_id
""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('--indir', metavar='PATH', type=str, required=True,
		help="""path to input snps directory for one species (contains files 'snps_*.txt')
requires having run 'merge_midas.py snps'""")
	parser.add_argument('--out', metavar='PATH', type=str, required=True,
		help="""path to output file containing list of markers""")
	parser.add_argument('--sample_map', metavar='PATH', type=str,
		help="""path to mapping file that specifes groups of related sample_ids.
the file should be tab-delimited with no header and two columns
example:
   sample1 subject1
   sample2 subject1
   sample3 subject2
   sample4 subject2
   sample5 subject3
""")
	parser.add_argument('--min_freq', type=float, metavar='FLOAT', default=0.10,
		help="""minimum allele frequency (proportion of reads) per site for SNP calling (0.10)""")
	parser.add_argument('--min_reads', type=int, metavar='INT', default=3,
		help="""minimum number of reads supporting allele per site for SNP calling (3)""")
	parser.add_argument('--max_groups', type=int, metavar='INT', default=1,
		help="""discard alleles found in > MAX_GROUPS (1)
groups are specified by the second field in '--sample_map'
if '--sample_map' is not specified, each sample_id is treated as a separate group""")
	parser.add_argument('--max_sites', type=int, metavar='INT',
		help="""maximum number of genomic sites to process (use all)
useful for quick tests""")
	args = vars(parser.parse_args())
	if args['sample_map']: args['samples'] = dict([_.rstrip().split() for _ in open(args['sample_map'])])
	else: args['samples'] = dict([(_.split()[0],_.split()[0]) for _ in open('%s/snps_summary.txt' % args['indir'])])
	return args

def track_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: track rare SNPs between all pairs of samples and determine transmission

Usage: strain_tracking.py track_markers [options]
""",
	epilog="""
Examples:
1) Track marker alleles for species_id 57955
strain_tracking.py track_markers --indir merged_snps/57955 --markers 57955.markers --out 57955.marker_sharing

2) Run a quick test for species_id 57955
strain_tracking.py track_markers --indir merged_snps/57955 --markers 57955.markers --out 57955.marker_sharing --max_sites 1000

Output fields:
  sample1: identifier for sample 1
  sample2: identifier for sample 2
  count1: number of marker alleles found in sample 1
  count2: number of marker alleles found in sample 2
  count_both: number of marker alleles found in sample 1 and 2
  count_either: number of marker alleles found in sample 1 or 2
""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('--indir', metavar='PATH', type=str, required=True,
		help="""path to input snps directory for one species (contains files 'snps_*.txt')
requires having run 'merge_midas.py snps'""")
	parser.add_argument('--out', metavar='PATH', type=str,
		help="""path to output file with marker sharing between all sample-pairs""")
	parser.add_argument('--markers', metavar='PATH', type=str,
		help="""path to list of marker alleles output by 'strain_tracking.py id_markers'""")
	parser.add_argument('--min_freq', type=float, metavar='FLOAT', default=0.10,
		help="""minimum allele frequency (proportion of reads) per site for SNP calling (0.10)""")
	parser.add_argument('--min_reads', type=int, metavar='INT', default=3,
		help="""minimum number of reads supporting allele per site for SNP calling (3)""")
	parser.add_argument('--max_sites', type=int, metavar='INT',
		help="""maximum number of sites to process (use all)
useful for quick tests""")
	parser.add_argument('--max_samples', type=int, metavar='INT',
		help="""maximum number of samples to process (use all)
useful for quick tests""")

	args = vars(parser.parse_args())
	if not os.path.isdir(args['indir']): sys.exit("Specified input directory '%s' does not exist" % args['indir'])
	if not os.path.isfile(args['markers']): sys.exit("Specified input file '%s' does not exist" % args['markers'])
	return args

def run_program(program, args):
	""" Run program specified by user (species, genes, or snps) """
	from midas.analyze import track_strains
	if program == 'id_markers':
		track_strains.id_markers(args)
	elif program == 'track_markers':
		track_strains.track_markers(args)
	else:
		sys.exit("Unrecognized program: '%s'" % program)

if __name__ == '__main__':
	program = get_program()
	args = get_arguments(program)
	utility.print_copyright()
	run_program(program, args)
