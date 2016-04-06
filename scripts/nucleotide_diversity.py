#!/usr/bin/env python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

# Next: allow user to specify subset of samples in input
#

import argparse, sys, os, gzip, numpy as np
from phylo_cnv import utility


class Sample:
	""" Base class for samples """
	def __init__(self, id):
		self.id = id
		self.total_pi = 0
		self.sites = 0
		self.snps = 0
		self.pi = 'NA'
		self.pnps = 'NA'
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
		self.total_pi += pi(freq)
		self.site_type[info['site_type']] += 1
		if is_snp(freq, maf): self.snp_type[info['site_type']] += 1

	def compute_diversity(self):
		self.sites = sum(self.site_type.values())
		self.snps = sum(self.snp_type.values())
		self.pi = self.total_pi/self.sites if self.sites > 0 else 'NA'
		self.pnps = self.compute_pnps()
		for gene in self.genes.values():
			gene.compute_diversity()

	def compute_pnps(self):
		if self.site_type['1D'] == 0 or self.snp_type['4D'] == 0 or self.site_type['4D'] == 0:
			return 'NA'
		else:
			pn = float(self.snp_type['1D'])/self.site_type['1D']
			ps = float(self.snp_type['4D'])/self.site_type['4D']
			return pn/ps

	def format_record(self, fields):
		values = {}
		values['sample_id'] = self.id
		values['pi'] = self.pi
		values['pnps'] = self.pnps
		values['sites'] = self.sites
		values['snps'] = self.snps
		for site_type in ['NC','1D','2D','3D','4D']:
			values['sites_'+site_type] = self.site_type[site_type]
		for snp_type in ['NC','1D','2D','3D','4D']:
			values['snps_'+snp_type] = self.snp_type[snp_type]
		record = '\t'.join([str(values[field]) for field in fields])+'\n'
		return record

class Gene:
	""" Base class for genes """
	def __init__(self, id):
		self.id = id
		self.total_pi = 0
		self.sites = 0
		self.snps = 0
		self.pi = 'NA'
		self.pnps = 'NA'
		self.site_type = {'NC':0,'1D':0,'2D':0,'3D':0,'4D':0}
		self.snp_type = {'NC':0,'1D':0,'2D':0,'3D':0,'4D':0}

	def update_stats(self, info, freq, maf):
		self.total_pi += pi(freq)
		self.site_type[info['site_type']] += 1
		if is_snp(freq, maf): self.snp_type[info['site_type']] += 1

	def compute_diversity(self):
		self.sites = sum(self.site_type.values())
		self.snps = sum(self.snp_type.values())
		self.pi = self.total_pi/self.sites if self.sites > 0 else 'NA'
		self.pnps = self.compute_pnps()
		
	def compute_pnps(self):
		if self.site_type['1D'] == 0 or self.snp_type['4D'] == 0 or self.site_type['4D'] == 0:
			return 'NA'
		else:
			pn = float(self.snp_type['1D'])/self.site_type['1D']
			ps = float(self.snp_type['4D'])/self.site_type['4D']
			return pn/ps

	def format_record(self, sample_id, fields):
		values = {}
		values['sample_id'] = sample_id
		values['gene_id'] = self.id
		values['pi'] = self.pi
		values['pnps'] = self.pnps
		values['sites'] = self.sites
		values['snps'] = self.snps
		for site_type in ['NC','1D','2D','3D','4D']:
			values['sites_'+site_type] = self.site_type[site_type]
		for snp_type in ['NC','1D','2D','3D','4D']:
			values['snps_'+snp_type] = self.snp_type[snp_type]
		record = '\t'.join([str(values[field]) for field in fields])+'\n'
		return record

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: nucleotide_diversity.py outdir [options]

Description: 
	use core-genome SNPs to estimate population genetic parameters of a species within and between samples
	within-sample parameters are estimated using the frequency of SNPs within individual samples
	between-sample parameters are estimated using the the mean frequency of SNPs across samples
Input: 
	SNP allele frequency matrix (SPECIES.snps.ref_freq), SNP annotation file (SPECIES.snps.info)
Output: 
	genome-wide nucleotide diversity, within and between samples
	per-gene nucleotide diversity, within and between samples
	genome-wide non-synonymous to synonymous SNP ratio (pN/pS), within and between samples
	per-gene non-synonymous to synonymous SNP ratio (pN/pS), within and between samples
""",
		epilog="""Examples:
1) Run using defaults:
nucleotide_diversity.py /path/to/outdir --ref_freq /path/to/SPECIES.snps.ref.freq --snp_info /path/to/SPECIES.snps.info

2) Use alternate definition of a SNP as 1% minor allele frequency:
nucleotide_diversity.py /path/to/outdir --maf 0.01 --ref_freq /path/to/SPECIES.snps.ref.freq --snp_info /path/to/SPECIES.snps.info

3) run a quick test:
nucleotide_diversity.py /path/to/outdir --max_sites 10000 --ref_freq /path/to/SPECIES.snps.ref.freq --snp_info /path/to/SPECIES.snps.info

""")
	parser.add_argument('outdir', type=str,
		help="""path to output directory""")
	parser.add_argument('--ref_freq', metavar='PATH', type=str, required=True,
		help="""input to reference frequency matrix""")
	parser.add_argument('--snp_info', metavar='PATH', type=str, required=True,
		help="""input to snp info file""")
	parser.add_argument('--max_sites', metavar='INT', type=int, required=False,
		help="""maximum number of sites to use from input matrix
useful for quick tests (use all)""")
	parser.add_argument('--maf', type=float, metavar='FLOAT', required=False, default=0.05,
		help="""minor allele frequency for SNP definition (0.05)""")

	args = vars(parser.parse_args())
	check_args(args)
	return args

def check_args(args):
	""" Check validity of command line arguments """
	pass

def parse_matrix(inpath):
	""" Parse SNP matrix """
	ext = inpath.split('.')[-1]
	infile = gzip.open(inpath) if ext == 'gz' else open(inpath)
	samples = next(infile).rstrip().split()[2:]
	for line in infile:
		values = line.rstrip().split()
		snp = values[0:2]
		freqs = values[2:]
		yield snp, samples, freqs

def fetch_samples(inpath):
	""" Get list of samples from SNP matrix """
	ext = inpath.split('.')[-1]
	infile = gzip.open(inpath) if ext == 'gz' else open(inpath)
	samples = next(infile).rstrip().split()[2:]
	return(samples)

def pi(x):
	""" Nucleotide diversity """
	return(2*x*(1-x))

def init_samples(args):
	""" Initialize samples """
	samples = {}
	for sample_id in fetch_samples(args['ref_freq']):
		samples[sample_id] = Sample(sample_id)
	samples['between'] = Sample('between')
	return samples

def is_snp(freq, min_maf):
	""" Determine if freq classified as a SNP or not """
	if min(freq, 1-freq) >= min_maf:
		return True
	else:
		return False

def parse_snps(args):
	""" Parse reference frequency matrix and snp info file """
	index = 0
	ref_freq = parse_matrix(args['ref_freq'])
	snp_info = utility.parse_file(args['snp_info'])
	while True:
		try:
			snp, sample_ids, freqs = next(ref_freq)
			info = next(snp_info)
		except StopIteration:
			break
		if [info['ref_id'], info['ref_pos']] != snp:
			sys.exit("Error: reference frequency matrix and SNP info file contain different sites")
		if args['max_sites'] and index >= args['max_sites']:
			break
		else:
			index += 1
			yield info, snp, sample_ids, freqs

def compute_diversity(args, samples):
	""" Compute number of sites and snps """
	# count sites and snps
	for info, snp, sample_ids, freqs in parse_snps(args):
		# within-sample
		for sample_id, freq in zip(sample_ids, freqs):
			if freq == 'NA': continue
			freq = float(freq)
			sample = samples[sample_id]
			sample.init_gene(info)
			sample.update_stats(info, freq, args['maf'])
		# between-sample
		if len(freqs) > 1:
			freq = np.mean([float(freq) for freq in freqs if freq != 'NA'])
			sample = samples['between']
			sample.init_gene(info)
			sample.update_stats(info, freq, args['maf'])
	# compute diversity metrics
	for sample in samples.values():
		sample.compute_diversity()

def write_per_sample_results(args, samples):
	""" Write per-sample results to output file(s) """
	outpath = ('%s/per-sample.diversity' % args['outdir'])
	outfile = open(outpath, 'w')
	fields = ['sample_id', 'pi', 'pnps', 'sites', 'snps',
			  'sites_NC', 'sites_1D', 'sites_2D', 'sites_3D', 'sites_4D',
			  'snps_NC', 'snps_1D', 'snps_2D', 'snps_3D', 'snps_4D']
	outfile.write('\t'.join(fields)+'\n')
	for sample in samples.values():
		record = sample.format_record(fields)
		outfile.write(record)

def write_per_gene_results(args, samples):
	""" Write per-gene results to output file(s) """
	outpath = ('%s/per-gene.diversity' % args['outdir'])
	outfile = open(outpath, 'w')
	fields = ['sample_id', 'gene_id', 'pi', 'pnps', 'sites', 'snps',
			  'sites_NC', 'sites_1D', 'sites_2D', 'sites_3D', 'sites_4D',
			  'snps_NC', 'snps_1D', 'snps_2D', 'snps_3D', 'snps_4D']
	outfile.write('\t'.join(fields)+'\n')
	gene_ids = set([])
	for sample in samples.values():
		for gene_id in sample.genes:
			gene_ids.add(gene_id)
	for sample in samples.values():
		for gene_id in gene_ids:
			if gene_id not in sample.genes:
				gene = Gene(gene_id)
			else:
				gene = sample.genes[gene_id]
			record = gene.format_record(sample.id, fields)
			outfile.write(record)

if __name__ == '__main__':
	args = parse_arguments()
	utility.print_copyright()
	print("Loading samples...\n")
	samples = init_samples(args)
	print("Estimating diversity metrics...\n")
	compute_diversity(args, samples)

	print("Writing results...\n")
	write_per_sample_results(args, samples)
	write_per_gene_results(args, samples)

