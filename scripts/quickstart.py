#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.2'

import argparse, sys, os, run_phylo_cnv

def parse_arguments():
	""" Get arguments for metagenomic species profiling """
	parser = argparse.ArgumentParser(usage='run_phylo_cnv.py species [options]')
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('-v', '--verbose', action='store_true', default=False)
	parser.add_argument('-i', type=str, dest='file', help='File containing list of paths to input metagenomes', required=True)
	parser.add_argument('-t', dest='threads', default=1, help='Number of threads to use')
	args = vars(parser.parse_args())
	return args


for program in ['species', 'genes', 'snvs']:
	args = get_arguments(program)
	check_arguments(program, args)
	run_program(program, args)



#def species_arguments():
#	""" Get arguments for metagenomic species profiling """
#	parser = argparse.ArgumentParser(usage='run_phylo_cnv.py species [options]')
#	parser.add_argument('program', help=argparse.SUPPRESS)
#	parser.add_argument('-v', '--verbose', action='store_true', default=False)
#	parser.add_argument('-1', type=str, dest='m1', help='FASTA/FASTQ file containing 1st mate if paired or unpaired reads', required=True)
#	parser.add_argument('-2', type=str, dest='m2', help='FASTA/FASTQ file containing 2nd mate if paired')
#	parser.add_argument('-o', type=str, dest='out', help='Path to output file', required=True)
#	parser.add_argument('-k', dest='keep_temp', default=False, action='store_true', help='Keep temporary files, including BLAST output')
#	parser.add_argument('-m', action='store_true', default=False, dest='norm', help='Estimate cellular relative abundance. Requires running MicrobeCensus and takes 20-30 minutes longer to complete.')
#	parser.add_argument('-s', type=str, dest='speed', default='fast', choices=['fast','sensitive'], help='Alignment speed/sensitivity (fast)')
#	parser.add_argument('-n', type=int, dest='reads', help='# reads to use from input file(s) (use all)')
#	parser.add_argument('-t', dest='threads', default=1, help='Number of threads to use')
#	args = vars(parser.parse_args())
#	return args




