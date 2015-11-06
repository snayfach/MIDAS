#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.2'

import os, sys
import argparse
import numpy as np

def print_copyright():
	print ("")
	print ("PhyloCNV: species abundance and strain-level genomic variation from metagenomes")
	print ("version %s; github.com/snayfach/PhyloCNV" % __version__)
	print ("Copyright (C) 2015 Stephen Nayfach")
	print ("Freely distributed under the GNU General Public License (GPLv3)")
	print ("")

def parse_arguments():
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__),
			description="""Merge species abundance files across samples. 
				Outputs include: a relative abundance matrix, a genome-coverage matrix, 
				and a table summarizing species prevalence and abundance across samples""")
	parser.add_argument('-v', '--verbose', action='store_true', default=False, help='verbose')
	parser.add_argument('-i', dest='indir', type=str, required=True,
		help="""input directory of species abundance files. filenames should have the format: <indir>/<sample_id>.species""")
	parser.add_argument('-o', dest='outbase', type=str, required=True,
		help="""basename for output files: <outbase>.species_abundance, <outbase>.species_coverage, <outbase>.species_prevalence""")
	parser.add_argument('-m', dest='min_cov', type=float, required=False, default=1.0,
		help="""minimum genome-coverage for estimating species prevalence (1.0)""")
	args = vars(parser.parse_args())
	add_annotations(args)
	check_args(args)
	return args

def print_arguments(args):
	print ("-------------------------------------------------------")
	print ("Merge Species Parameters:")
	print ("Input directory: %s" % args['indir'])
	print ("Output basename: %s" % args['outbase'])
	print ("Minimum coverage for estimating prevalence: %s" % args['min_cov'])
	print ("-------------------------------------------------------")
	print ("")

def check_args(args):
	if not os.path.isdir(args['indir']):
		sys.exit("Input directory (-i) not found:\n%s" % os.path.abspath(args['indir']))
	base_dir = os.path.dirname(os.path.abspath(args['outbase']))
	if not os.path.isdir(base_dir):
		sys.exit("Output directory of base name (-o) not found:\n%s" % base_dir)
	if not os.path.isfile(args['annotations']):
		sys.exit("Species annotation file not found:\n%s" % os.path.abspath(args['annotations']))

def add_annotations(args):
	maindir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	args['annotations'] = '%s/ref_db/annotations.txt' % maindir

def parse_phylo_species(inpath):
	infile = open(inpath)
	names = next(infile).rstrip().split()
	for line in infile:
		values = line.rstrip().split()
		yield dict([(name,value) for name,value in zip(names,values)])

def list_samples(args):
	samples = []
	for file in os.listdir(args['indir']):
		x = file.split('.')
		if len(x) > 1 and x[-1] == 'species':
			samples.append('.'.join(x[0:-1]))
	if len(samples) == 0:
		sys.exit("No <sample_id>.species files found in:\n%s" % os.path.abspath(indir))
	return(samples)

def list_species(args, samples):
	species_ids = []
	inpath = '%s/%s.species' % (args['indir'], samples[0])
	for r in parse_phylo_species(inpath):
		species_ids.append(r['cluster_id'])
	return species_ids

def store_data(args, sample_ids, species_ids):
	data = {}
	for species in species_ids:
		data[species] = {'abundance':[], 'coverage':[]}
	for sample in sample_ids:
		inpath = '%s/%s.species' % (args['indir'], sample)
		for r in parse_phylo_species(inpath):
			data[r['cluster_id']]['abundance'].append(float(r['relative_abundance']))
			data[r['cluster_id']]['coverage'].append(float(r['coverage']))
	return data

def prevalence(x, y):
	return(sum([1 if _ >= y else 0 for _ in x]))

def compute_stats(args, data):
	species_ids = data.keys()
	stats = dict([(s,{}) for s in species_ids])
	for species in species_ids:
		# abundance
		x = data[species]['abundance']
		stats[species]['median_abundance'] = np.median(x)
		stats[species]['mean_abundance'] = np.mean(x)
		# coverage
		x = data[species]['coverage']
		stats[species]['median_coverage'] = np.median(x)
		stats[species]['mean_coverage'] = np.mean(x)
		# prevalence
		stats[species]['prevalence'] = prevalence(x, y=args['min_cov'])
	return stats

def read_annotations(args):
	annotations = {}
	infile = open(args['annotations'])
	fields = next(infile).rstrip().split('\t')
	for line in infile:
		values = line.rstrip().split('\t')
		record = dict([(f,v) for f,v in zip(fields, values)])
		annotations[record['cluster_id']] = record['consensus_name']
	return annotations

def write_abundance(args, sample_ids, species_ids, data):
	outfile = open('%s.species_abundance' % args['outbase'], 'w')
	outfile.write('\t'.join(['species_id']+sample_ids)+'\n')
	for species in species_ids:
		outfile.write(species)
		for x in data[species]['abundance']:
			outfile.write('\t%s' % str(x))
		outfile.write('\n')

def write_coverage(args, sample_ids, species_ids, data):
	outfile = open('%s.species_coverage' % args['outbase'], 'w')
	outfile.write('\t'.join(['species_id']+sample_ids)+'\n')
	for species in species_ids:
		outfile.write(species)
		for x in data[species]['coverage']:
			outfile.write('\t%s' % str(x))
		outfile.write('\n')

def write_stats(args, sample_ids, species_ids, stats, annotations):
	# open output file
	outfile = open('%s.species_prevalence' % args['outbase'], 'w')
	fields = ['mean_coverage', 'median_coverage', 'mean_abundance', 'median_abundance', 'prevalence']
	header = ['species_id', 'species_name']+fields
	outfile.write('\t'.join(header)+'\n')
	# order species ids by decreasing abundance
	sorted_species = [[x, y['prevalence']] for x, y in stats.items()]
	sorted_species.sort(key=lambda x: x[1], reverse=True)
	# write stats
	for species, prevalence in sorted_species:
		outfile.write(species)
		outfile.write('\t%s' % annotations[species])
		for field in fields:
			outfile.write('\t%s' % str(round(stats[species][field], 2)))
		outfile.write('\n')

if __name__ == "__main__":

	args = parse_arguments()
	if args['verbose']: print_copyright()
	if args['verbose']: print_arguments(args)
	
	# list samples and species
	sample_ids = list_samples(args)
	species_ids = list_species(args, sample_ids)
	annotations = read_annotations(args)
	
	# read in data & compute stats
	data = store_data(args, sample_ids, species_ids)
	stats = compute_stats(args, data)

	# write results
	write_abundance(args, sample_ids, species_ids, data)
	write_coverage(args, sample_ids, species_ids, data)
	write_stats(args, sample_ids, species_ids, stats, annotations)




