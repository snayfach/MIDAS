#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, sys
from midas import utility

class Species:
	""" Base class for species """
	def __init__(self, id, info):
		self.id = id
		self.consensus_name = info[id]['consensus_name']
		self.samples = []

class Sample:
	""" Base class for samples """
	def __init__(self, dir):
		self.dir = dir
		self.id = os.path.basename(self.dir)
		self.paths = self.init_paths()

	def init_paths(self):
		paths = {}
		species = '/'.join([self.dir, 'species/species_profile.txt'])
		if os.path.isfile(species): paths['species'] = species
		else: paths['species'] = None
		snps = '/'.join([self.dir, 'snps/summary.txt'])
		if os.path.isfile(snps): paths['snps'] = snps
		else: paths['snps'] = None
		genes = '/'.join([self.dir, 'genes/summary.txt'])
		if os.path.isfile(genes): paths['genes'] = genes
		else: paths['genes'] = None
		return paths

def write_summary_stats(species_id, samples, args, type):
	""" Write summary file for samples """
	outfile = open('%s/%s/%s.%s.summary' % (args['outdir'], species_id, species_id, type), 'w')
	if type == 'snps':
		fields = ['genome_length', 'covered_bases', 'fraction_covered', 'mean_coverage']
	else:
		fields = ['pangenome_size', 'covered_genes', 'fraction_covered', 'mean_coverage', 'marker_coverage']
	outfile.write('\t'.join(['sample_id']+fields)+'\n')
	for sample in samples:
		stats = read_stats(sample.paths[type], type)
		outfile.write(sample.id)
		for field in fields:
			if field in stats[species_id]:
				outfile.write('\t' + str(stats[species_id][field]))
			else: # temporary fix
				outfile.write('\t')
		outfile.write('\n')

def sort_species(species):
	""" Sort list of species by number of samples in descending order """
	x = sorted([(sp, len(sp.samples)) for sp in species], key=lambda x: x[1], reverse=True)
	return [_[0] for _ in x]

def select_species(args, type='genes'):
	""" Select all species with a minimum number of high-coverage samples"""
	# read species annotations
	species_info = {}
	inpath = os.path.join(args['db'], 'species_info.txt')
	for rec in utility.parse_file(inpath):
		species_info[rec['species_id']] = rec
	# fetch all species with at least 1 sample
	species = {}
	for sample in load_samples(args):
		if not sample.paths[type]:
			sys.stderr.write("Warning: no genes output for sample: %s\n" % sample.dir)
			continue
		for id, info in read_stats(sample.paths[type], type).items():
			if (args['species_id']
					and id not in args['species_id'].split(',')):
				continue # skip unspecified species
			elif (args['max_samples']
					and id in species
					and len(species[id].samples) >= args['max_samples']):
				continue # skip species with too many samples
			elif float(info['mean_coverage']) < args['sample_depth']:
				continue # skip low-coverage sample
			elif type=='snps' and float(info['fraction_covered']) < args['fract_cov']:
				continue # skip low-coverage sample
			if id not in species:
				species[id] = Species(id, species_info) # initialize new species
			species[id].samples.append(sample) # append sample
	# dict to list
	species = species.values()
	# remove species with an insufficient number of samples
	species = [sp for sp in species if len(sp.samples) >= args['min_samples']]
	# sort by number of samples
	species = sort_species(species)
	# select a subset of species to analyze
	if args['max_species'] is not None and len(species) > args['max_species']:
		species = species[0:args['max_species']]
	# create species directories
	for sp in species:
		outdir = os.path.join(args['outdir'], sp.id)
		if not os.path.isdir(outdir): os.mkdir(outdir)
	print "  found %s species with sufficient high-coverage samples\n" % len(species)
	return species

def read_stats(inpath, type):
	stats = {}
	for rec in utility.parse_file(inpath):
		if 'cluster_id' in rec: rec['species_id'] = rec['cluster_id']
		if 'phyeco_coverage' in rec: rec['marker_coverage'] = rec['phyeco_coverage']
		if 'average_depth' in rec: rec['mean_coverage'] = rec['average_depth']
		if type=='genes':
			rec['fraction_covered'] = float(rec['covered_genes'])/float(rec['pangenome_size'])
		stats[rec['species_id']] = rec
	return stats

def list_samples(input, intype):
	if intype == 'dir':
		return([os.path.join(input, _) for _ in os.listdir(input)])
	elif intype == 'file':
		return([x.rstrip() for x in open(input).readlines()])
	elif intype == 'list':
		return(input.split(','))

def load_samples(args):
	samples = []
	for dir in list_samples(args['input'], args['intype']):
		if os.path.isdir(dir):
			samples.append(Sample(dir))
		else:
			sys.stderr.write("Warning: specified sample dir does not exist: %s\n" % dir)
	return samples

