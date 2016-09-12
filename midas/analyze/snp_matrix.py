#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, numpy as np, random
from midas import utility

class GenomicSite:
	""" Base class for genomic sites """
	def __init__(self, files, samples, info=None):
		try:
			self.ref_freq = next(files['ref_freq'])[1]
			self.depth = next(files['depth'])[1]
			self.alt_allele = next(files['alt_allele'])[1]
			self.info = next(info)
			self.id = self.info['site_id']
			self.ref_id, self.ref_pos, self.ref_allele = self.parse_id()
			self.samples = samples
			self.cast_numeric()
			self.hq_samples = len([_ for _ in samples if _.pass_qc])
			
		except StopIteration:
			self.id = None

	def filter_site(self, site_prev, site_maf, site_types):
		""" determine if site passes quality control """
		if self.prev < max(1e-6, site_prev): # must be data for at least 1 sample
			return True
		elif self.maf < site_maf:
			return True
		elif self.ref_allele not in ['A','T','C','G']:
			return True
		elif self.info['site_type'] not in site_types:
			return True
		else:
			return False

	def summary_stats(self, weighted):
		self.prev = float(len(self.samples))/self.hq_samples
		mean_freq = self.mean_freq(weighted)
		self.maf = min(mean_freq, 1-mean_freq)

	def mean_freq(self, weighted):
		if weighted:
			sum_depth = sum(site.depth)
			weighted_freqs = [freq * depth / sum_depth for freq, depth in zip(self.ref_freq, self.depth)]
			return np.mean(weighted_freqs)
		else:
			return np.mean(self.ref_freq)
		
	def pooled_pi(self, site_depth, weighted=False):
		x = self.maf(site_depth, weighted)
		if x == 'NA': return 0
		else: return 2*x*(1-x)

	def prune(self, indexes):
		""" prune low quality samples """
		self.ref_freq = [self.ref_freq[i] for i in indexes]
		self.depth = [self.depth[i] for i in indexes]
		self.alt_allele = [self.alt_allele[i] for i in indexes]
		self.samples = [self.samples[i] for i in indexes]
				
	def prune_samples(self, site_depth, site_ratio):
		""" remove data for samples with low coverage at site """
		index = 0
		indexes = []
		for depth, sample in zip(self.depth, self.samples):
			if (sample.pass_qc and
					depth >= site_depth and
					depth/sample.depth <= site_ratio):
				indexes.append(index)
			index += 1

		self.prune(indexes)

	def resample_reads(self, rand_reads):
		index = 0
		for ref_freq, depth in zip(self.ref_freq, self.depth):
			if depth > rand_reads:
				self.depth[index] = rand_reads
				if ref_freq > 0 and ref_freq < 1:
					count_ref = int(round(ref_freq * depth))
					count_alt = depth - count_ref
					alleles = random.sample([1]*count_ref+[0]*count_alt, rand_reads)
					self.ref_freq[index] = np.mean(alleles)
			index += 1

	def resample_samples(self, rand_samples):
		indexes = random.sample(range(len(self.samples)), rand_samples)
		self.ref_freq = [self.ref_freq[i] for i in indexes]
		self.depth = [self.depth[i] for i in indexes]
		self.samples = [self.samples[i] for i in indexes]
		
	def cast_numeric(self):
		self.ref_freq = [float(_) if _ != 'NA' else 0.0 for _ in self.ref_freq ]
		self.depth = [int(_) if _ != 'NA' else 0.0 for _ in self.depth]
	
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

def parse_sites(indir, samples):
	""" yield genomic sites from input files """
	index = 0
	files = {} # open input files
	for ext in ['alt_allele', 'depth', 'ref_freq']:
		files[ext] = parse_tsv('%s/snps_%s.txt' % (indir, ext))
	info = utility.parse_file('%s/snps_info.txt' % indir)

	while True: # yield GenomicSite
		site = GenomicSite(files, samples, info)
		if not site.id:
			break
		else:
			index += 1
			yield site
	for file in files.values(): # close input files
		file.close()

class Sample:
	""" Base class for sample """
	def __init__(self, info, filters, stop):
		self.id = info['sample_id']
		self.depth = float(info['mean_coverage'])
		self.pass_qc = self.qc(info, filters, stop)

	def qc(self, info, filters, stop):
		if stop:
			return False
		elif float(info['fraction_covered']) < filters['min_cov']:
			return False
		elif float(info['mean_coverage']) < filters['min_depth']:
			return False
		else:
			return True

def fetch_samples(args):
	samples = []
	count_hq = 0
	stop = False
	filters = {'min_depth':args['sample_depth'], 'min_cov':args['fract_cov']}
	for info in utility.parse_file('%s/snps_summary.txt' % args['indir']):
		samples.append(Sample(info, filters, stop))
		if samples[-1].pass_qc: count_hq += 1
		if count_hq >= args['max_samples']: stop = True
	return samples

