#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3

import argparse, sys
from phylo_cnv import phylo_species

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(description='Estimate the abundance of genome-clusters from metagenomes.')
	
	io = parser.add_argument_group('Input/Output (required)')
	io.add_argument('-i', type=str, dest='inpath', required=True,
		help='path to input metagenome in FASTQ/FASTA format. gzip (.gz) and bzip (.bz2) compression supported')
	io.add_argument('-o', type=str, dest='outpath', required=True,
		help='output file name')
	
	speed = parser.add_argument_group('Pipeline Speed (optional)')
	speed.add_argument('-n', type=int, dest='nreads', required=False, default=float('Inf'),
		help='number of reads to use from input metagenome (use all)')
	speed.add_argument('-p', type=int, dest='threads', required=False, default=1,
		help='number of threads to use for database search (1)')
	
	qc = parser.add_argument_group('Quality control (optional)')
	qc.add_argument('-q', type=int, dest='min_quality', required=False, default=0,
		help='keep reads with quality >= MIN_QUALITY (0)')
	qc.add_argument('-l', type=int, dest='min_length', required=False, default=0,
		help='keep reads with length >= MIN_LENGTH (0)')
	qc.add_argument('-u', type=float, dest='max_n', required=False, default=1.0,
		help='keep reads with fraction unknown bases <= MAX_N (1.0)')
	
	parser.add_argument('-v', dest='verbose', action='store_true', default=False)

	return parser.parse_args()

if __name__ == '__main__':
	args = vars(parse_arguments())
	cluster_abundance = phylo_species.estimate_species_abundance(args)
	phylo_species.write_abundance(args['outpath'], cluster_abundance)
