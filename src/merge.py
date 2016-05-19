# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, random, utility, sys

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


def select_species(args, type='genes'):
	""" Select all species with a minimum number of high-coverage samples"""
	# read species annotations
	species_info = {}
	inpath = os.path.join(args['db'], 'annotations.txt')
	for rec in utility.parse_file(inpath):
		if 'cluster_id' in rec: rec['species_id'] = rec['cluster_id'] ## temporary fix
		species_info[rec['species_id']] = rec
	# fetch all species with at least 1 sample
	species = {}
	for sample in load_samples(args):
		if not sample.paths[type]:
			sys.stderr.write("Warning: no genes output for sample: %s\n" % sample.dir)
			continue
		for id, info in read_stats(sample.paths[type]).items():
			if (args['species_id']
					and id not in args['species_id'].split(',')):
				continue # skip unspecified species
			elif (args['max_samples']
					and id in species
					and len(species[id].samples) >= args['max_samples']):
				continue # skip species with too many samples
			elif (type == 'genes'
					and float(info['mean_coverage']) < args['sample_depth']):
				continue # skip low-coverage species
			elif (type == 'snps'
					and float(info['average_depth']) < args['sample_depth']):
				continue # skip low-coverage species
			if id not in species:
				species[id] = Species(id, species_info) # initialize new species
			species[id].samples.append(sample) # append sample
	# remove species with an insufficient number of samples
	for sp in species.copy():
		if len(species[sp].samples) < args['min_samples']:
			del species[sp]
	species = species.values()
	# select a subset of species to analyze
	if args['max_species'] is not None and len(species) > args['max_species']:
		species = random.sample(species, args['max_species'])
	print "  found %s species with sufficient high-coverage samples\n" % len(species)
	return species

def read_stats(inpath):
	stats = {}
	infile = open(inpath)
	fields = next(infile).rstrip().split('\t')[1:]
	for line in infile:
		species_id = line.rstrip().split('\t')[0]
		values = line.rstrip().split('\t')[1:]
		stats[species_id] = dict([(i,j) for i,j in zip(fields, values)])
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

