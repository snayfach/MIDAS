#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, numpy as np
from midas import utility

class GenomicSite:
	""" Base class for genomic sites """
	def __init__(self, files, samples, info=None):
		try:
			self.id, self.ref_freq = next(files['ref_freq'])
			self.id, self.depth = next(files['depth'])
			self.id, self.alt_allele = next(files['alt_allele'])
			self.info = next(info) if info else None

			self.ref_id, self.ref_pos, self.ref_allele = self.parse_id()
			self.samples = samples
			
		except StopIteration:
			self.id = None

	def parse_id(self):
		ref_id = self.id.rsplit('|', 2)[0]
		ref_pos = int(self.id.split('|')[2])
		ref_allele = self.id.split('|')[3] if len(self.id.split('|')) >= 4 else '' # testing
		return ref_id, ref_pos, ref_allele

	def sample_values(self): # return a dic mapping sample id to site values
		d = {}
		for index, sample in enumerate(self.samples):
			d[sample] = {}
			d[sample]['ref_freq'] = self.ref_freq[index]
			d[sample]['depth'] = self.depth[index]
			d[sample]['alt_allele'] = self.alt_allele[index]
		return d

	def prev(self, site_depth):
		count = len([_ for _ in self.depth if int(_) >= site_depth])
		return float(count)/len(self.depth)

	def mean_freq(self):
		freqs = [float(_) for _ in self.ref_freq if _ != 'NA']
		if len(freqs) > 0: return np.mean(freqs)
		else: return 'NA'

	def mean_depth(self):
		depths = [float(_) for _ in self.depth if _ != 'NA']
		if len(depths) > 0: return np.mean(depths)
		else: return 'NA'

	def maf(self):
		x = self.mean_freq()
		if x == 'NA': return x
		else: return min(x, 1-x)

	def allele_props(self):
		sums = {'A':0.0, 'T':0.0, 'C':0.0, 'G':0.0}
		props = {}
		# sum alternate allele proportions
		for freq, allele in zip(self.ref_freq, self.alt_allele):
			if allele != 'NA': sums[allele] += 1-float(freq)
		# sum reference allele proportions
		for freq in self.ref_freq:
			if freq != 'NA': sums[self.ref_allele] += float(freq)
		# normalize proportions
		total = sum(sums.values())
		if total > 0:
			for allele in sums:
				props[allele] = float('%.3g' % (sums[allele]/total))
		return props

	def filter(self, site_depth, site_prev, site_maf):
		if self.prev(site_depth) < site_prev:
			return True
		elif self.maf() < site_maf:
			return True
		elif self.ref_allele not in ['A','T','C','G']:
			return True
		else:
			return False

def init_paths(indir):
	""" fetch paths to input files """
	paths = {}
	exts = ['alt_allele', 'info', 'summary', 'depth', 'ref_freq']
	for ext in exts:
		inpath = '%s/snps_%s.txt' % (indir, ext)
		if os.path.isfile(inpath):
			paths[ext] = inpath
	return paths

def parse_tsv(inpath):
	""" yield records from tab-delimited file with header """
	infile = utility.iopen(inpath)
	header = next(infile).rstrip('\n').split('\t')
	for line in infile:
		split_line = line.rstrip('\n').split('\t')
		id = split_line[0]
		values = split_line[1:]
		yield id, values
	infile.close()

def parse_sites(indir):
	""" yield genomic sites from input files """
	index = 0
	files = {} # open input files
	for ext, path in init_paths(indir).items():
		files[ext] = parse_tsv(path)
	samples = list_samples(indir)
	while True: # yield GenomicSite
		site = GenomicSite(files, samples)
		if not site.id:
			break
		else:
			index += 1
			yield site
	for file in files.values(): # close input files
		file.close()

def check_sample_ids(args, paths):
	""" make sure that specified sample ids are present in input """
	samples = [_.split()[0] for _ in open(paths['summary'])][1:]
	for sample_id in args['samples'].copy():
		if sample_id not in samples:
			del args['samples'][sample_id]
	if len(args['samples']) == 0:
		sys.exit("No samples remain")

def list_samples(indir, max_samples=None):
	""" list sample ids from specified input """
	infile = open('%s/snps_ref_freq.txt' % indir)
	sample_ids = next(infile).rstrip('\n').split('\t')[1:]
	if max_samples is not None: sample_ids = sample_ids[0:max_samples]
	return sample_ids

