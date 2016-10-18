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

	def cast_numeric(self):
		""" format numeric fields """
		self.ref_freq = [float(_) if _ != 'NA' else 0.0 for _ in self.ref_freq ]
		self.depth = [int(_) if _ != 'NA' else 0.0 for _ in self.depth]
	
	def parse_id(self):
		""" parse info from site_id """
		ref_id = self.id.rsplit('|')[0]
		ref_pos = int(self.id.split('|')[1])
		ref_allele = self.id.split('|')[2]
		return ref_id, ref_pos, ref_allele
		
	def filter_site(self, site_prev, site_maf, site_type, site_freq):
		""" determine if site passes quality control """
		if self.prev < max(1e-6, site_prev):
			return True, 'site_prev'
		elif self.maf < site_maf:
			return True, 'site_maf'
		elif self.ref_allele not in ['A','T','C','G']:
			return True, 'ref_allele'
		elif self.info['site_type'] not in site_type:
			return True, 'site_type'
		elif self.refalt_freq < site_freq:
			return True, 'site_freq'
		else:
			return False, None

	def summary_stats(self, weighted):
		""" compute summary stats for site across hq samples """
		self.prev = float(len(self.samples))/self.hq_samples
		self.mean_freq = self.compute_mean_freq(weighted)
		self.maf = min(self.mean_freq, 1-self.mean_freq)
		self.allele_freqs = self.compute_allele_freqs()
		self.refalt_freq = self.compute_refalt_freq()
		
	def compute_allele_freqs(self):
		""" compute combined frequency of 4 alleles across samples; only considers 2 alleles per sample """
		freqs = {'A':0.0, 'T':0.0, 'C':0.0, 'G':0.0}
		for freq, alt_allele in zip(self.ref_freq, self.alt_allele):
			freqs[self.ref_allele] += float(freq)
			if alt_allele != 'NA': freqs[alt_allele] += 1-float(freq)
		total = sum(freqs.values())
		if total > 0:
			for allele in freqs:
				freqs[allele] = freqs[allele]/total
		return freqs
		
	def compute_refalt_freq(self):
		""" compute combined frequency of reference and major alternate allele """
		freqs = self.allele_freqs
		ref_freq = freqs[self.ref_allele]
		alt_freq = 0
		for allele, freq in self.allele_freqs.items():
			if allele != self.ref_allele:
				alt_freq = max(freq, alt_freq)
		return ref_freq + alt_freq

	def compute_mean_freq(self, weighted):
		""" compute average frequency of reference allele with optional weighting of samples """
		if len(self.ref_freq) == 0:
			return 0.0
		elif weighted:
			sum_depth = sum(site.depth)
			weighted_freqs = [freq * depth / sum_depth for freq, depth in zip(self.ref_freq, self.depth)]
			return np.mean(weighted_freqs)
		else:
			return np.mean(self.ref_freq)

	def prune(self, indexes):
		""" subset site based on list of indexes """
		self.ref_freq = [self.ref_freq[i] for i in indexes]
		self.depth = [self.depth[i] for i in indexes]
		self.alt_allele = [self.alt_allele[i] for i in indexes]
		self.samples = [self.samples[i] for i in indexes]
				
	def prune_samples(self, site_depth, site_ratio):
		""" filter samples at site based on coverage """
		index = 0
		indexes = []
		for depth, sample in zip(self.depth, self.samples):
			if (sample.pass_qc and
					depth >= site_depth and
					depth/sample.depth <= site_ratio):
				indexes.append(index)
			index += 1
		self.prune(indexes)

	def resample_reads(self, rand_reads, replace_reads):
		""" resample random number of reads per sample """
		index = 0
		for ref_freq, depth in zip(self.ref_freq, self.depth):
			self.depth[index] = rand_reads
			if ref_freq > 0 and ref_freq < 1:
				count_ref = int(round(ref_freq * depth))
				count_alt = depth - count_ref
				alleles = np.random.choice([1]*count_ref+[0]*count_alt, rand_reads, replace=replace_reads)
				self.ref_freq[index] = np.mean(alleles)
			index += 1

def parse_tsv(inpath):
	""" yield records from tab-delimited file with row and column names """
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
	def __init__(self, info):
		self.id = info['sample_id']
		self.info = info
		self.depth = float(self.info['mean_coverage'])
		self.covered = float(self.info['fraction_covered'])

	def filter(self, sample_depth, fract_cov):
		if self.covered < fract_cov:
			return True
		elif self.depth < sample_depth:
			return True
		else:
			return False

def fetch_samples(args):
	""" initialize samples from indir and set pass_qc flag """
	samples = []
	count_hq = 0
	stop = False
	filters = {'min_depth':args['sample_depth'], 'min_cov':args['fract_cov']}
	for info in utility.parse_file('%s/snps_summary.txt' % args['indir']):
		samples.append(Sample(info, filters, stop))
		if samples[-1].pass_qc: count_hq += 1
		if count_hq >= args['max_samples']: stop = True
	return samples

