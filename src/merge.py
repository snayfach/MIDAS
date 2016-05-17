# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os

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
		snps = '/'.join([self.dir, 'snps/snps_summary_stats.txt'])
		if os.path.isfile(snps): paths['snps'] = snps
		else: paths['snps'] = None
		genes = '/'.join([self.dir, 'genes/genes_summary_stats.txt'])
		if os.path.isfile(genes): paths['genes'] = genes
		else: paths['genes'] = None
		return paths

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
		if not os.path.isdir(input):
			sys.exit("\nSpecified input directory does not exist:\n%s" % input)
		else:
			return([os.path.join(input, _) for _ in os.listdir(input)])
	elif intype == 'file':
		if not os.path.isfile(input):
			sys.exit("\nSpecified input file does not exist:\n%s" % input)
		else:
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

