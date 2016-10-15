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
Usage: build_midas_db.py [options]
""")
	parser.add_argument('--type', choices=['marker', 'gene', 'genome', 'all'], default='all',
		help="""Type of database to build (all)
marker: build database of marker genes for metagenomic species profiling
        enables you to use the command "run_midas.py species"
  gene: build database of genes for metagenomic pan-genome profiling
        enables you to use the command "run_midas.py genes"
genome: build database of representative genomes for metagenomic SNP profiling
        enables you to use the command "run_midas.py snps"
   all: build all three types of databases sequentially""")
	parser.add_argument('--outdir', type=str, required=True,
		help="Path to directory to store database files")
	parser.add_argument('--genomes', type=str, required=True,
		help="""Path to directory of input genomes
each subdirectory should contain files for one genome and should be named according to the genome_id""")
	parser.add_argument('--mapping', type=str, required=True,
		help="""Path to tab-delimited genome mapping file with header
field 1: genome_id
         character field corresponding to subdirectory within --genomes directory
field 2: species_id
         character field indicating species identifier for genome_id
field 3: rep_genome
         optional boolean field (0/1) indicating if a genome_id is a representative for species_id""")
	parser.add_argument('--threads', type=str, default=1,
		help="Number of threads to use")
	parser.add_argument('--pid', type=float, default=0.95,
		help="Percent identity threshold for defining non-redundant genes (0.95)")
	parser.add_argument('--iter_size', type=int, default=500,
		help="""Maximum number of genomes to process with at one time
to prevent exceeding USEARCH's 4G memory limit (500)""")
	parser.add_argument('--max_species', type=int, default=float('inf'),
		help="Maximum number of species to process from input (use all).\nUseful for quick tests")
	parser.add_argument('--max_genomes', type=int, default=float('inf'),
		help="Maximum number of genomes to process per species (use all).\nUseful for quick tests")
	parser.add_argument('--compress', action='store_true', default=False,
		help="Compress output files with gzip")

	args = vars(parser.parse_args())
	
	return args

def check_args(args):

	if not os.path.isdir(args['outdir']):
		os.mkdir(args['outdir'])
	if not os.path.isdir(args['genomes']):
		sys.exit("\nError: could not locate directory specified by --genomes: %s" % args['genomes'])
	if not os.path.isfile(args['mapping']):
		sys.exit("\nError: could not locate file specified by --mapping: %s" % args['mapping'])

	for program in ['hmmsearch', 'usearch']:
		if not utility.which(program):
			error = ""
			error += "\nError: program '%s' not found in your PATH" % program
			error += "\nMake sure that you've installed the program and added it's location to your PATH"
			sys.exit(error)


if __name__ == "__main__":

	args = fetch_arguments()
	utility.add_executables(args)
	check_args(args)
	utility.print_copyright()
	build_db.run_pipeline(args)






