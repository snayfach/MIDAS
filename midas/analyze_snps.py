
# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, sys, numpy as np
from midas import utility

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

	def init_gene(self, gene_id):
		if gene_id not in self.genes:
			self.genes[gene_id] = Gene(gene_id)
	
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
	def __init__(self, sample_ids, freq, depth, alleles, site_depth=0, info=None):
		try:
			rec_freq = next(freq)
			self.id = rec_freq.split('\t')[0]
			self.ref_id = self.id.rsplit('|', 1)[0]
			self.ref_allele = self.id.split('|')[-1]
			self.ref_pos = int(self.id.split('|')[-2])
			self.info = next(info) if info else None
			self.freqs = rec_freq.rstrip().split('\t')[1:]
			self.depths = next(depth).rstrip().split('\t')[1:]
			self.alleles = next(alleles).rstrip().split('\t')[1:]
			self.sample_ids = sample_ids
			self.count_samples = len(self.sample_ids)
			self.mean_depth = self.mean_depth()
			self.mean_freq = self.mean_freq()
			self.prev = self.prev(site_depth)
			self.maf = maf(self.mean_freq)
		except StopIteration:
			self.id = None

	def mean_freq(self):
		freqs = [float(_) for _ in self.freqs if _ != 'NA']
		if len(freqs) > 0: return np.mean(freqs)
		else: return 'NA'
		
	def mean_depth(self):
		depths = [int(_) for _ in self.depths if _ != 'NA']
		if len(depths) > 0: return np.mean(depths)
		else: return 'NA'

	def prev(self, site_depth):
		count = len([_ for _ in self.depths if int(_) >= site_depth])
		return float(count)/self.count_samples

	def allele_props(self):
		sums = {'A':0.0, 'T':0.0, 'C':0.0, 'G':0.0}
		props = {}
		# sum alternate allele proportions
		for freq, allele in zip(self.freqs, self.alleles):
			if allele != 'NA': sums[allele] += 1-float(freq)
		# sum reference allele proportions
		for freq in self.freqs:
			if freq != 'NA': sums[self.ref_allele] += float(freq)
		# normalize proportions
		total = sum(sums.values())
		if total > 0:
			for allele in sums:
				props[allele] = sums[allele]/total
		return props

	def iterate(self, sample_ids):
		for sample_id, freq, count in zip(self.sample_ids, self.freqs, self.depths):
			if sample_id not in sample_ids:
				continue
			elif freq == 'NA':
				continue
			else:
				yield sample_id, float(freq), int(count)

def parse_sites(indir, site_depth=0, info_path=None, max_sites=None):
	""" Parse reference frequency matrix and snp info file """
	index = 0
	sample_ids = list_sample_ids('%s/snps_ref_freq.txt' % indir)
	freq = open('%s/snps_ref_freq.txt' % indir); next(freq)
	depth = open('%s/snps_depth.txt' % indir); next(depth)
	alleles = open('%s/snps_alt_allele.txt' % indir); next(alleles)
	info = utility.parse_file(info_path) if info_path else None
	while True:
		site = GenomicSite(sample_ids, freq, depth, alleles, site_depth, info)
		if not site.id:
			break
		elif max_sites and index >= max_sites:
			break
		else:
			index += 1
			yield site
	freq.close(); depth.close()
	
def init_paths(args):
	paths = {}
	for ext in ['ref_freq', 'alt_allele', 'depth', 'info', 'summary']:
		inpath = '%s.snps.%s' % (args['inbase'], ext)
		if os.path.isfile(inpath):
			paths[ext] = inpath
		else:
			sys.exit('Input file does not exist: %s' % inpath)
	return paths

def add_summary_stats(args, samples):
	""" add summary stats """
	for rec in utility.parse_file(args['summary']):
		if 'run_accession' in rec: rec['sample_id'] = rec['run_accession']
		if rec['sample_id'] in samples:
			samples[rec['sample_id']].mean_depth = float(rec['average_depth'])
			samples[rec['sample_id']].fraction_covered = float(rec['fraction_covered'])
	samples['between'].mean_depth = np.mean([_.mean_depth for _ in samples.values() if _.id != 'between'])

def list_sample_ids(freq_path):
	""" Fetch list of site, sample ids from file """
	infile = open(freq_path)
	ids = next(infile).rstrip().split('\t')
	infile.close()
	return ids[1:]

def init_samples(args):
	""" Initialize samples """
	samples = {}
	for sample_id in list_sample_ids(args):
		if not args['samples'] or sample_id in args['samples'].split(','):
			samples[sample_id] = Sample(sample_id)
	samples['between'] = Sample('between')
	add_summary_stats(args, samples)
	return samples

def maf(x):
	""" Minor allele frequency """
	if x == 'NA':
		return x
	else:
		return min(x, 1-x)

def is_snp(freq, min_maf):
	""" Determine if freq classified as a SNP or not """
	if maf(freq) >= min_maf:
		return True
	else:
		return False

def pi(x):
	""" Nucleotide diversity """
	return(2*x*(1-x))

def identify_snps(args, samples): ## Re-write
	""" Store snps per sample """
	for site in parse_sites(args):
		freqs = []
		for sample_id, freq, count in site.iterate(samples):
			freqs.append(freq)
			if is_snp(freq, args['maf']):
				samples[sample_id].snps.add(site.id)
		if len(freqs) > 1 and is_snp(np.mean(freqs), args['maf']):
			samples['between'].snps.add(site.id)
