#!/usr/bin/env python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip, numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from phylo_cnv import utility

class Sample:
	""" Base class for samples """
	def __init__(self, id):
		self.id = id
		self.snps, self.sites, self.pi_bp, self.pnps_pi, self.pnps_snps = 5*['NA']
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
		self.sites = sum(self.site_type.values())
		self.snps = sum(self.snp_type.values())
		self.pi_bp = sum(self.pi_type.values())/self.sites if self.sites > 0 else 'NA'
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
		values['sites'] = self.sites
		values['snps'] = self.snps
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
		self.snps, self.sites, self.pi_bp, self.pnps_pi, self.pnps_snps = 5*['NA']
		self.pi_type = {'NC':0,'1D':0,'2D':0,'3D':0,'4D':0}
		self.site_type = {'NC':0,'1D':0,'2D':0,'3D':0,'4D':0}
		self.snp_type = {'NC':0,'1D':0,'2D':0,'3D':0,'4D':0}

	def update_stats(self, info, freq, maf):
		self.site_type[info['site_type']] += 1
		self.pi_type[info['site_type']] += pi(freq)
		if is_snp(freq, maf): self.snp_type[info['site_type']] += 1

	def compute_diversity(self):
		self.sites = sum(self.site_type.values())
		self.snps = sum(self.snp_type.values())
		self.pi_bp = sum(self.pi_type.values())/self.sites if self.sites > 0 else 'NA'
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
		values['sites'] = self.sites
		values['snps'] = self.snps
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
	parser.add_argument('--inbase', metavar='PATH', type=str, required=True,
		help="""basename for input files""")
	parser.add_argument('--max_sites', metavar='INT', type=int, required=False,
		help="""maximum number of sites to use from input matrix
useful for quick tests (use all)""")
	parser.add_argument('--maf', type=float, metavar='FLOAT', required=False, default=0.05,
		help="""minor allele frequency for SNP definition (0.05)""")
	parser.add_argument('--samples', type=str, metavar='STR', required=False,
		help="""comma-separated list of samples to use for computing diversity metrics (use all)""")
	args = vars(parser.parse_args())
	args.update(init_paths(args))
	return args

def init_paths(args):
	paths = {}
	for ext in ['ref_freq', 'alt_allele', 'depth', 'info', 'sample_ids', 'summary']:
		inpath = '%s.snps.%s' % (args['inbase'], ext)
		if os.path.isfile(inpath):
			paths[ext] = inpath
		else:
			sys.exit('Input file does not exist: %s' % inpath)
	return paths

def init_samples(args):
	""" Initialize samples """
	samples = {}
	for line in open(args['sample_ids']):
		sample_id = line.rstrip()
		if not args['samples'] or sample_id in args['samples'].split(','):
			samples[sample_id] = Sample(sample_id)
	samples['between'] = Sample('between')
	return samples

def is_snp(freq, min_maf):
	""" Determine if freq classified as a SNP or not """
	if min(freq, 1-freq) >= min_maf:
		return True
	else:
		return False

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

def pi(x):
	""" Nucleotide diversity """
	return(2*x*(1-x))

def compute_diversity(args, samples):
	""" Compute number of sites and snps """
	# count sites and snps
	for site in parse_sites(args):
		freqs_filt = []
		# within-sample
		for sample_id, freq in zip(site.samples, site.freqs):
			if freq != 'NA' and sample_id in samples:
				samples[sample_id].init_gene(site.info)
				samples[sample_id].update_stats(site.info, float(freq), args['maf'])
				freqs_filt.append(float(freq))
		# between-sample
		if len(freqs_filt) > 1:
			samples['between'].init_gene(site.info)
			samples['between'].update_stats(site.info, np.mean(freqs_filt), args['maf'])
	# compute diversity metrics
	for sample in samples.values():
		sample.compute_diversity()

def write_per_sample_results(args, samples):
	""" Write per-sample results to output file(s) """
	outpath = ('%s/per-sample.diversity' % args['outdir'])
	outfile = open(outpath, 'w')
	fields = ['sample_id', 'pi_bp', 'pnps_snps', 'pnps_pi', 'sites', 'snps',
			  'sites_NC', 'sites_1D', 'sites_2D', 'sites_3D', 'sites_4D',
			  'snps_NC', 'snps_1D', 'snps_2D', 'snps_3D', 'snps_4D',
			  'pi_NC', 'pi_1D', 'pi_2D', 'pi_3D', 'pi_4D']
	outfile.write('\t'.join(fields)+'\n')
	for sample in samples.values():
		record = sample.format_record(fields)
		outfile.write(record)

def write_per_gene_results(args, samples):
	""" Write per-gene results to output file(s) """
	outpath = ('%s/per-gene.diversity' % args['outdir'])
	outfile = open(outpath, 'w')
	fields = ['sample_id', 'gene_id', 'pi_bp', 'pnps_snps', 'pnps_pi', 'sites', 'snps',
			  'sites_NC', 'sites_1D', 'sites_2D', 'sites_3D', 'sites_4D',
			  'snps_NC', 'snps_1D', 'snps_2D', 'snps_3D', 'snps_4D',
			  'pi_NC', 'pi_1D', 'pi_2D', 'pi_3D', 'pi_4D']
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
	
	print("Loading samples...")
	samples = init_samples(args)
	
	print("Estimating diversity metrics...")
	compute_diversity(args, samples)

	print("Writing results...")
	write_per_sample_results(args, samples)
	write_per_gene_results(args, samples)

