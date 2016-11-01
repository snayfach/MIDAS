#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, gzip, numpy as np, random, csv
from midas.utility import print_copyright

class Sample:
	""" Base class for sample """
	def __init__(self, info):
		self.id = info['sample_id']
		self.info = info
		self.mean_depth = float(self.info['mean_coverage'])
		self.fract_cov = float(self.info['fraction_covered'])

	def filter(self, mean_depth, fract_cov):
		if self.fract_cov < fract_cov:
			return True
		elif self.mean_depth < mean_depth:
			return True
		else:
			return False

class Species:
	""" Base class for species """
	
	def __init__(self, dir):
		
		self.dir = dir
		self.id = os.path.basename(self.dir)
		self.init_paths()
		self.init_files()
		self.init_samples()

	def init_paths(self):
		self.paths = {}
		for type in ['alt_allele', 'depth', 'info', 'ref_freq', 'summary']:
			self.paths[type] = '%s/snps_%s.txt' % (self.dir, type)

	def init_files(self):
		self.files = {}
		for type in ['alt_allele', 'depth', 'info', 'ref_freq', 'summary']:
			file = open(self.paths[type])
			if type in ['info', 'summary']:
				self.files[type] = csv.DictReader(file, delimiter='\t')
			else:
				self.files[type] = csv.reader(file, delimiter='\t')

	def init_samples(self):
		self.sample_ids = None
		for file in ['alt_allele', 'depth', 'ref_freq']:
			self.sample_ids = next(self.files[file])[1:]


class GenomicSite:
	""" Base class for genomic sites """
	def __init__(self, species, samples):
		try:
			self.samples = samples
			self.fetch_row(species)
			self.info = next(species.files['info'])
			self.id=1
			self.site_id()
		except StopIteration:
			self.id = None

	def site_id(self):
		""" parse info from site_id """
		self.id = self.info['site_id']
		self.ref_id, self.ref_pos, self.ref_allele = self.id.rsplit('|', 2)
		
	def fetch_row(self, species):
		""" fetch next row from matrices """
		self.sample_ids = species.sample_ids
		values = zip(self.sample_ids,
		             next(species.files['ref_freq'])[1:],
					 next(species.files['depth'])[1:],
					 next(species.files['alt_allele'])[1:]
					 )
		for sample_id, ref_freq, depth, alt_allele in values:
			try:
				sample = self.samples[sample_id]
				sample.ref_freq = float(ref_freq) if ref_freq != 'NA' else 0.0
				sample.depth = int(depth)
				sample.alt_allele = alt_allele
			except KeyError:
				pass
				
	def filter_samples(self, site_depth, site_ratio):
		""" filter samples at site based on coverage """
		for sample in self.samples.values():
			if (sample.depth >= site_depth and
					sample.depth/sample.mean_depth <= site_ratio):
				sample.keep = True
			else:
				sample.keep = False

	def filter_site(self, site_prev=None, site_maf=None, site_type=None, site_freq=None):
		""" determine if site passes quality control """
		if self.ref_allele not in ['A','T','C','G']:
			self.keep = False
		if site_prev and self.prev < max(1e-6, site_prev):
			self.keep = False
		elif site_maf and self.minor_freq < site_maf:
			self.keep = False
		elif site_type and self.info['site_type'] not in site_type:
			self.keep = False
		else:
			self.keep = True

	def prevalence(self):
		self.count_samples = sum([1 for s in self.samples.values() if s.keep])
		return float(self.count_samples)/len(self.samples)

	def summary_stats(self, weighted=False):
		""" compute summary stats for site across samples """
		self.prev = self.prevalence()
		self.mean_freq = self.compute_mean_freq(weighted)
		self.minor_freq = min(self.mean_freq, 1-self.mean_freq)

	def compute_mean_freq(self, weighted):
		""" compute average frequency of reference allele with optional weighting of samples """
		if self.count_samples == 0:
			return 0.0
		else:
			return np.mean([s.ref_freq for s in self.samples.values() if s.keep])

def fetch_samples(species, mean_depth=0, fract_cov=0, max_samples=float('inf'), keep_samples=None, exclude_samples=None):
	""" List and select samples from input """
	samples = {}
	for info in species.files['summary']:
		# init sample
		sample = Sample(info)
		# excludes sample
		if sample.filter(mean_depth, fract_cov):
			continue
		if keep_samples and sample.id not in keep_samples:
			continue
		if exclude_samples and sample.id in exclude_samples:
			continue
		if len(samples) >= max_samples:
			continue
		# store sample
		samples[sample.id] = sample
	if len(samples) == 0:
		sys.exit("\nError: no samples satisfied your selection criteria.\nTry running again with more lenient parameters\n")
	return samples


def fetch_sites(species, samples):
	""" yield genomic sites from species across samples """
	index = 0
	while True:
		site = GenomicSite(species, samples)
		if not site.id:
			break
		else:
			index += 1
			yield site