#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, sys, csv
from midas import utility

class Species:
	""" Base class for species """
	def __init__(self, id, species_info, genome_info):
		self.id = id
		self.samples = []
		self.info = species_info[self.id]
		self.genome_info = genome_info[self.info['rep_genome']]

	def fetch_sample_depth(self):
		self.sample_depth = []
		for sample in self.samples:
			self.sample_depth.append(float(sample.info[self.id]['mean_coverage']))

	def write_sample_info(self, dtype, outdir):
		""" Write summary file for samples """
		outfile = open('%s/%s/%s_summary.txt' % (outdir, self.id, dtype), 'w')
		if dtype == 'snps':
			fields = ['genome_length', 'covered_bases', 'fraction_covered', 'mean_coverage']
		else:
			fields = ['pangenome_size', 'covered_genes', 'fraction_covered', 'mean_coverage', 'marker_coverage']
		outfile.write('\t'.join(['sample_id']+fields)+'\n')
		for sample in self.samples:
			path = '%s/%s/summary.txt' % (sample.dir, dtype)
			outfile.write(sample.id)
			for field in fields:
				value = sample.info[self.id][field]
				outfile.write('\t' + str(value))
			outfile.write('\n')

	def open_outfiles(self, dtype, outdir):
		""" Open output files for species """
		self.outdir = os.path.join(outdir, self.id)
		self.files = {}
		if dtype == 'snps':
			for ftype in ['info', 'freq', 'depth']:
				path = '%s/snps_%s.txt' % (self.outdir, ftype)
				self.files[ftype] = open(path, 'w')
			for type in ['freq', 'depth']:
				record = ['site_id']+[s.id for s in self.samples]
				self.files[ftype].write('\t'.join(record)+'\n')
			info_fields = ['site_id', 'site_prev', 'ref_allele', 'major_allele',
						   'minor_allele', 'minor_freq', 'atcg_counts', 'site_type',
						   'atcg_aas', 'gene_id']
			self.files['info'].write('\t'.join(info_fields)+'\n')

	def close_outfiles(self):
		for file in self.files.values():
			file.close()

class Sample:
	""" Base class for samples """
	def __init__(self, dir, data_type):
		self.dir = dir
		self.id = os.path.basename(self.dir)
		self.info = self.read_info(data_type)

	def read_info(self, data_type):
		path = '%s/%s/summary.txt' % (self.dir, data_type)
		if os.path.isfile(path):
			info = {}
			for r in csv.DictReader(open(path), delimiter='\t'):
				info[r['species_id']] = r
			return info
		else:
			return None

def init_samples(indirs, data_type):
	""" Initialize samples """
	samples = []
	for dir in indirs:
		sample = Sample(dir, data_type)
		if sample.info is None:
			pass
			#sys.stderr.write("Warning: missing/incomplete output: %s\n" % dir)
		else:
			samples.append(sample)
	return samples

def read_species_info(db):
	""" Read species annotations """
	species_info = {}
	path = os.path.join(db, 'species_info.txt')
	for r in csv.DictReader(open(path), delimiter='\t'):
		species_info[r['species_id']] = r
	return species_info

def read_genome_info(db):
	""" Read genome annotations """
	genome_info = {}
	path = os.path.join(db, 'genome_info.txt')
	for r in csv.DictReader(open(path), delimiter='\t'):
		genome_info[r['genome_id']] = r
	return genome_info

def filter_sample_species(sample, species, species_id, args, dtype):
	""" Determine whether sample-species pair fails filters """
	info = sample.info[species_id]
	if (args['species_id']
			and species_id not in args['species_id'].split(',')):
		return True # skip unspecified species
	elif (args['max_samples']
			and species_id in species
			and len(species[species_id].samples) >= args['max_samples']):
		return True # skip species with too many samples
	elif float(info['mean_coverage']) < args['sample_depth']:
		return True # skip low-coverage sample
	elif dtype=='snps' and float(info['fraction_covered']) < args['fract_cov']:
		return True # skip low-coverage sample
	else:
		return False

def sort_species(species):
	""" Sort list of species by number of samples in descending order """
	x = sorted([(sp, len(sp.samples)) for sp in species], key=lambda x: x[1], reverse=True)
	return [_[0] for _ in x]

def init_species(samples, args, dtype):
	""" Store high quality sample-species pairs """
	species = {}
	species_info = read_species_info(args['db'])
	genome_info = read_genome_info(args['db'])
	for sample in samples:
		for species_id in sample.info:
			if species_id not in species:
				species[species_id] = Species(species_id, species_info, genome_info)
			if filter_sample_species(sample, species, species_id, args, dtype):
				continue
			else:
				species[species_id].samples.append(sample)
	return species.values()

def filter_species(species, args):
	""" Pick subset of species to analyze """
	keep = []
	for sp in sort_species(species):
		sp.nsamples = len(sp.samples)
		if sp.nsamples < int(args['min_samples']):
			continue
		elif args['max_species'] and len(keep) >= args['max_species']:
			continue
		else:
			sp.fetch_sample_depth()
			sp.outdir = args['outdir']+'/'+sp.id
			keep.append(sp)
			if not os.path.isdir(sp.outdir):
				os.mkdir(sp.outdir)
	return keep

def select_species(args, dtype):
	""" Select all species with a minimum number of high-coverage samples"""
	samples = init_samples(args['indirs'], dtype)
	species = init_species(samples, args, dtype)
	species = filter_species(species, args)
	return species

def write_snps_readme(args, sp):
	outfile = open('%s/%s/readme.txt' % (args['outdir'], sp.id), 'w')
	outfile.write("""
Description of output files and file formats from 'merge_midas.py snps'

Output files
############
snps_freq.txt
  frequency of minor allele per genomic site and per sample
  see: snps_info.txt for major, minor, and reference alleles
snps_depth.txt
  number of reads mapped to genomic site per sample
  only accounts for reads matching either major or minor allele
snps_info.txt  
  metadata for genomic site
  see below for more information
snps_summary.txt
  alignment summary statistics per sample
  see below for more information
snps_log.txt
  log file containing parameters used

Output formats
############
snps_freq.txt and snps_depth.txt
  tab-delimited matrix files
  field names are sample ids
  row names are genome site ids
  genomic sites have the format: ref_id|ref_pos (eg: NC_ATTCG|1 corresponds to position 1 on contig NC_ATTCG)
snps_summary.txt
  sample_id: sample identifier
  genome_length: number of base pairs in representative genome
  covered_bases: number of reference sites with at least 1 mapped read
  fraction_covered: proportion of reference sites with at least 1 mapped read
  mean_coverage: average read-depth across reference sites with at least 1 mapped read
snps_info.txt
  site_id: genomic site_id, format=ref_id|ref_pos
  ref_id: scaffold/contig identifier
  ref_pos: position of site on ref_id
  ref_allele: allele in reference genome
  major_allele: most common allele in metagenomes
  minor_allele: second most common allele in metagenomes
  count_samples: number of metagenomes where site_id was found
  count_atcg: counts of all 4 alleles (A,T,C,G) in pooled metagenomes
  snp_type: site is [MONO,BI,TRI,QUAD]-allelic
  site_type: NC (non-coding), 1D, 2D, 3D, 4D (degeneracy)
  amino_acid_atcg: amino acids encoded by 4 possible alleles
  gene_id: gene that intersects site

Additional information for species can be found in the reference database:
 %s/rep_genomes/%s
""" % (args['db'], sp.id) )
	outfile.close()

