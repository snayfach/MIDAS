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
		self.genome_info = genome_info[self.info['representative_genome']]

	def fetch_sample_depth(self):
		self.sample_depth = []
		for sample in self.samples:
			self.sample_depth.append(float(sample.info[self.id]['mean_coverage']))

	def write_sample_info(self, dtype, outdir):
		""" Write summary file for samples """
		outfile = open('%s/%s/%s_summary.txt' % (outdir, self.id, dtype), 'w')
		if dtype == 'snps':
			fields = ['genome_length', 'covered_bases', 'fraction_covered', 'mean_coverage', 'aligned_reads', 'mapped_reads']
		else:
			fields = ['pangenome_size', 'covered_genes', 'fraction_covered', 'mean_coverage', 'marker_coverage', 'aligned_reads', 'mapped_reads']
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
			sys.stderr.write("Warning: missing/incomplete output: %s\n" % dir)
		else:
			samples.append(sample)
	return samples

def read_species_info(iggdb):
	""" Read species annotations """
	return iggdb.species

def read_genome_info(iggdb):
	""" Read genome annotations """
	return iggdb.genomes

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
	species_info = read_species_info(args['iggdb'])
	genome_info = read_genome_info(args['iggdb'])
	for sample in samples:
		for species_id in sample.info:
			if species_id not in species:
				species[species_id] = Species(species_id, species_info, genome_info)
			if not filter_sample_species(sample, species, species_id, args, dtype):
				species[species_id].samples.append(sample)
	return list(species.values())

def filter_species(species, args):
	""" Pick subset of species to analyze """
	keep = []
	for sp in sort_species(species):
		sp.nsamples = len(sp.samples)
		if sp.nsamples < int(args['min_samples']):
			continue
		if args['max_species'] and len(keep) >= args['max_species']:
			continue
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



