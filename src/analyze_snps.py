
# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, sys, utility, numpy as np

class Sample:
	""" Base class for samples """
	def __init__(self, id):
		self.id = id
		self.snps = set([])
		self.count_snps, self.count_sites = 0, 0
		self.pi_bp, self.pnps_pi, self.pnps_snps = 3*['NA']
		self.pi_type = {'NC':0,'1D':0,'2D':0,'3D':0,'4D':0}
		self.site_type = {'NC':0,'1D':0,'2D':0,'3D':0,'4D':0}
		self.snp_type = {'NC':0,'1D':0,'2D':0,'3D':0,'4D':0}
		self.genes = {}

	def init_gene(self, info):
		if info['gene_id'] == 'NA':
			return
		elif info['gene_id'] in self.genes:
			return
		else:
			self.genes[info['gene_id']] = Gene(info['gene_id'])
	
	def update_stats(self, info, freq, maf):
		if info['gene_id'] != 'NA': self.genes[info['gene_id']].update_stats(info, freq, maf)
		self.site_type[info['site_type']] += 1
		self.pi_type[info['site_type']] += pi(freq)
		if is_snp(freq, maf): self.snp_type[info['site_type']] += 1

	def compute_diversity(self):
		self.count_sites = sum(self.site_type.values())
		self.count_snps = sum(self.snp_type.values())
		self.pi_bp = sum(self.pi_type.values())/self.count_sites if self.count_sites > 0 else 'NA'
		self.pnps_snps = self.compute_pnps(type='snps')
		self.pnps_pi = self.compute_pnps(type='pi')
		for gene in self.genes.values():
			gene.compute_diversity()

	def compute_pnps(self, type, na_value='NA'):
		if self.site_type['1D'] == 0 or self.site_type['4D'] == 0:
			return na_value
		elif type == 'snps':
			if self.snp_type['4D'] == 0: return na_value
			pn = float(self.snp_type['1D'])/self.site_type['1D']
			ps = float(self.snp_type['4D'])/self.site_type['4D']
		elif type == 'pi':
			if self.pi_type['4D'] == 0: return na_value
			pn = float(self.pi_type['1D'])/self.site_type['1D']
			ps = float(self.pi_type['4D'])/self.site_type['4D']
		return pn/ps

	def format_record(self, fields):
		values = {}
		values['sample_id'] = self.id
		values['pi_bp'] = self.pi_bp
		values['pnps_snps'] = self.pnps_snps
		values['pnps_pi'] = self.pnps_pi
		values['sites'] = self.count_sites
		values['snps'] = self.count_snps
		for site_type in ['NC','1D','2D','3D','4D']:
			values['sites_'+site_type] = self.site_type[site_type]
		for snp_type in ['NC','1D','2D','3D','4D']:
			values['snps_'+snp_type] = self.snp_type[snp_type]
		for pi_type in ['NC','1D','2D','3D','4D']:
			values['pi_'+pi_type] = self.pi_type[pi_type]
		record = '\t'.join([str(values[field]) for field in fields])+'\n'
		return record

class Gene:
	""" Base class for genes """
	def __init__(self, id):
		self.id = id
		self.count_snps, self.count_sites = 0, 0
		self.pi_bp, self.pnps_pi, self.pnps_snps = 3*['NA']
		self.pi_type = {'NC':0,'1D':0,'2D':0,'3D':0,'4D':0}
		self.site_type = {'NC':0,'1D':0,'2D':0,'3D':0,'4D':0}
		self.snp_type = {'NC':0,'1D':0,'2D':0,'3D':0,'4D':0}

	def update_stats(self, info, freq, maf):
		self.site_type[info['site_type']] += 1
		self.pi_type[info['site_type']] += pi(freq)
		if is_snp(freq, maf): self.snp_type[info['site_type']] += 1

	def compute_diversity(self):
		self.count_sites = sum(self.site_type.values())
		self.count_snps = sum(self.snp_type.values())
		self.pi_bp = sum(self.pi_type.values())/self.count_sites if self.count_sites > 0 else 'NA'
		self.pnps_snps = self.compute_pnps(type='snps')
		self.pnps_pi = self.compute_pnps(type='pi')
		
	def compute_pnps(self, type, na_value='NA'):
		if self.site_type['1D'] == 0 or self.site_type['4D'] == 0:
			return na_value
		elif type == 'snps':
			if self.snp_type['4D'] == 0: return na_value
			pn = float(self.snp_type['1D'])/self.site_type['1D']
			ps = float(self.snp_type['4D'])/self.site_type['4D']
		elif type == 'pi':
			if self.pi_type['4D'] == 0: return na_value
			pn = float(self.pi_type['1D'])/self.site_type['1D']
			ps = float(self.pi_type['4D'])/self.site_type['4D']
		return pn/ps

	def format_record(self, sample_id, fields):
		values = {}
		values['sample_id'] = sample_id
		values['gene_id'] = self.id
		values['pi_bp'] = self.pi_bp
		values['pnps_snps'] = self.pnps_snps
		values['pnps_pi'] = self.pnps_pi
		values['sites'] = self.count_sites
		values['snps'] = self.count_snps
		for site_type in ['NC','1D','2D','3D','4D']:
			values['sites_'+site_type] = self.site_type[site_type]
		for snp_type in ['NC','1D','2D','3D','4D']:
			values['snps_'+snp_type] = self.snp_type[snp_type]
		for pi_type in ['NC','1D','2D','3D','4D']:
			values['pi_'+pi_type] = self.pi_type[pi_type]
		record = '\t'.join([str(values[field]) for field in fields])+'\n'
		return record

class GenomicSite:
	""" Base class for genomic sites """
	def __init__(self, info, samples, freq, depth):
		try:
			self.info = next(info)
			self.id = '|'.join([self.info['ref_id'], self.info['ref_pos']])
			self.freqs = next(freq).rstrip().split()[2:]
			self.depths = next(depth).rstrip().split()[2:]
			self.samples = samples
		except StopIteration:
			self.id = None

def init_paths(args):
	paths = {}
	for ext in ['ref_freq', 'alt_allele', 'depth', 'info', 'sample_ids', 'summary']:
		inpath = '%s.snps.%s' % (args['inbase'], ext)
		if os.path.isfile(inpath):
			paths[ext] = inpath
		else:
			sys.exit('Input file does not exist: %s' % inpath)
	return paths

def add_summary_stats(args, samples):
	""" add summary stats """
	for rec in utility.parse_file(args['summary']):
		id = rec['run_accession'] ## This field should be updated to sample_id in the future!!
		if id in samples:
			samples[id].mean_depth = float(rec['average_depth'])
	samples['between'].mean_depth = np.mean([_.mean_depth for _ in samples.values() if _.id != 'between'])

def init_samples(args):
	""" Initialize samples """
	samples = {}
	for line in open(args['sample_ids']):
		sample_id = line.rstrip()
		if not args['samples'] or sample_id in args['samples'].split(','):
			samples[sample_id] = Sample(sample_id)
	samples['between'] = Sample('between')
	add_summary_stats(args, samples)
	return samples

def is_snp(freq, min_maf):
	""" Determine if freq classified as a SNP or not """
	if freq == 'NA':
		return False
	freq = float(freq)
	if min(freq, 1-freq) >= min_maf:
		return True
	else:
		return False

def pi(x):
	""" Nucleotide diversity """
	return(2*x*(1-x))

def parse_sites(args):
	""" Parse reference frequency matrix and snp info file """
	index = 0
	samples = [_.rstrip() for _ in open(args['sample_ids'])]
	freq = open(args['ref_freq']); next(freq)
	depth = open(args['depth']); next(depth)
	info = utility.parse_file(args['info'])
	while True:
		site = GenomicSite(info, samples, freq, depth)
		if not site.id:
			break
		elif args['max_sites'] and index >= args['max_sites']:
			break
		else:
			index += 1
			yield site
	freq.close(); depth.close()

def identify_snps(args, samples):
	""" Store snps per sample """
	for site in parse_sites(args):
		for sample_id, freq, count in zip(site.samples, site.freqs, site.depths):
			if sample_id not in samples:
				continue
			if is_snp(freq, args['maf']):
				samples[sample_id].snps.add(site.id)
		freq = np.mean([float(_) for _ in site.freqs if _ != 'NA'])
		if is_snp(freq, args['maf']):
			samples['between'].snps.add(site.id)
