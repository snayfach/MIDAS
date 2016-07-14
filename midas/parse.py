#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os

class GenomicSite:
	""" Base class for genomic sites """
	def __init__(self, files):
		try:
			self.ref_freq = next(files['ref_freq'])
			self.depth = next(files['depth'])
			self.alt_allele = next(files['alt_allele'])
			self.info = next(files['info'])
			self.id = self.info['site_id']
			self.ref_id = self.id.rsplit('|', 1)[0]
			self.ref_pos = int(self.id.rsplit('|', 1)[1])
		except StopIteration:
			self.id = None

def parse_tsv(inpath):
	""" yield records from tab-delimited file with header """
	infile = open(inpath)
	fields = next(infile).rstrip().split('\t')
	for line in infile:
		values = line.rstrip().split('\t')
		if len(fields) == len(values):
			yield dict([(i,j) for i,j in zip(fields, values)])
	infile.close()

def parse_sites(args, paths):
	""" yield genomic sites from input files """
	index = 0
	files = {} # open input files
	for ext, path in paths.items():
		files[ext] = parse_tsv(path)
	while True: # yield GenomicSite
		site = GenomicSite(files)
		if not site.id:
			break
		elif args['max_sites'] is not None and index >= args['max_sites']:
			break
		else:
			index += 1
			yield site
	for file in files.values(): # close input files
		file.close()

