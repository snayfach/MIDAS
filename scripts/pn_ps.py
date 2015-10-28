#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)


# Recipe:
# Count number

__version__ = '0.0.2'

import argparse, sys, os, gzip, Bio.SeqIO

def parse_arguments():
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
	parser.add_argument('-i', dest='ref_freq', type=str, help='Input allele frequency matrix', required=True)
	return vars(parser.parse_args())





