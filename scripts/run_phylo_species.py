#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3

import argparse, sys
from phylo_cnv import phylo_species

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(description='Estimate the abundance of genome-clusters from metagenomes.')
	parser.add_argument('-i', type=str, dest='inpaths', required=True, nargs='+')
	parser.add_argument('-d', type=str, dest='indir', required=False)
	parser.add_argument('-o', type=str, dest='outbase', required=True)
	parser.add_argument('-n', type=int, dest='nreads', required=False, default=float('Inf'))
	parser.add_argument('-q', type=int, dest='min_quality', required=False, default=0)
	parser.add_argument('-l', type=int, dest='min_length', required=False, default=0)
	parser.add_argument('-u', type=float, dest='max_n', required=False, default=1.0)
	parser.add_argument('-t', type=str, dest='temp_dir', required=False)
	parser.add_argument('-v', dest='verbose', action='store_true', default=False)
	return parser.parse_args()

if __name__ == '__main__':
	args = vars(parse_arguments())
	cluster_abundance, cluster_summary = phylo_species.estimate_species_abundance(args)
	phylo_species.write_abundance(args['outbase']+'.abundance', cluster_abundance)
	phylo_species.write_summary(args['outbase']+'.summary', cluster_summary)
