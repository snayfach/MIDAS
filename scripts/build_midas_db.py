#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, subprocess, sys, argparse, shutil
from midas import utility
from midas.build import build_db
import Bio.SeqIO

def fetch_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description:
This script will allow you to build your own custom MIDAS database
Usage: build_midas_db.py indir mapfile outdir [options]
""")
	parser.add_argument('indir', type=str,
		help="""Path to directory of input genomes
Each subdirectory should be named according to a genome_id
Each subdirectory should contain the files:
  <genome_id>.fna: Genomic DNA sequence in FASTA format 
  <genome_id>.ffn: Gene DNA sequences in FASTA format 
  <genome_id>.features: Genomic coordinates of genes. Tab-delimited with a header fields:
     * scaffold_id (CHAR)
     * gene_id (CHAR)
     * start (INT)
     * end (INT)
     * strand ('+' or '-')	  
""")
	parser.add_argument('mapfile', type=str,
		help="""
Path to mapping file that specifies which genomes belonging to the same species.
The file should be tab-delimited file with a header and 3 fields:
  * genome_id (CHAR): corresponds to subdirectory within INDIR
  * species_id (CHAR): species identifier for genome_id
  * rep_genome (0 or 1): indicator if genome_id should be used for SNP calling
""")
	parser.add_argument('outdir', type=str,
		help="Directory to store MIDAS database")
	parser.add_argument('--threads', type=str, metavar='INT', default=1,
		help="Number of threads to use")
	parser.add_argument('--pid', type=float, metavar='FLOAT', default=0.95,
		help="Percent identity threshold for defining non-redundant genes (0.95)")
	parser.add_argument('--iter_size', type=int, default=500, metavar='INT',
		help="""Maximum number of genomes to process with at one time
to prevent exceeding USEARCH's 4G memory limit (500)""")
	parser.add_argument('--max_species', type=int, default=float('inf'), metavar='INT',
		help="Maximum number of species to process from input (use all).\nUseful for quick tests")
	parser.add_argument('--max_genomes', type=int, default=float('inf'), metavar='INT',
		help="Maximum number of genomes to process per species (use all).\nUseful for quick tests")
	parser.add_argument('--compress', action='store_true', default=False,
		help="Compress output files with gzip")

	args = vars(parser.parse_args())
	
	return args

def check_args(args):

	if not os.path.isdir(args['outdir']):
		os.mkdir(args['outdir'])
	if not os.path.isdir(args['indir']):
		sys.exit("\nError: could not locate directory specified by --genomes: %s\n" % args['indir'])
	if not os.path.isfile(args['mapfile']):
		sys.exit("\nError: could not locate file specified by --mapping: %s\n" % args['mapfile'])

	for program in ['hmmsearch', 'usearch']:
		if not utility.which(program):
			error = "\nError: program '%s' not found in your PATH" % program
			error += "\nMake sure that you've installed the program and added it's location to your PATH\n"
			sys.exit(error)


if __name__ == "__main__":
	args = fetch_arguments()
	utility.add_executables(args)
	check_args(args)
	utility.print_copyright()
	build_db.run_pipeline(args)






