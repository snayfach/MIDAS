#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, sys, numpy as np
from midas.run import species
from midas.merge import merge

def store_data(args, samples, species_ids):
	data = {}
	for species_id in species_ids:
		data[species_id] = {}
		for field in ['relative_abundance', 'coverage', 'count_reads']:
			data[species_id][field] = []
	for sample in samples:
		abundance = species.read_abundance(sample.paths['species'])
		for species_id, values in abundance.items():
			for field in ['relative_abundance', 'coverage', 'count_reads']:
				if field in values: # temporary fix
					data[species_id][field].append(values[field])
	return data

def prevalence(x, y):
	return(sum([1 if _ >= y else 0 for _ in x]))

def compute_stats(args, data):
	species_ids = data.keys()
	stats = dict([(s,{}) for s in species_ids])
	for species_id in species_ids:
		# abundance
		x = data[species_id]['relative_abundance']
		stats[species_id]['median_abundance'] = np.median(x)
		stats[species_id]['mean_abundance'] = np.mean(x)
		# coverage
		x = data[species_id]['coverage']
		stats[species_id]['median_coverage'] = np.median(x)
		stats[species_id]['mean_coverage'] = np.mean(x)
		# prevalence
		stats[species_id]['prevalence'] = prevalence(x, y=args['min_cov'])
	return stats

def write_abundance(args, samples, species_info, data):
	for field in ['relative_abundance', 'coverage', 'count_reads']:
		outfile = open('%s/%s.txt' % (args['outdir'], field), 'w')
		outfile.write('\t'.join(['species_id']+[s.id for s in samples])+'\n')
		for species_id in species_info:
			outfile.write(species_id)
			for x in data[species_id][field]:
				outfile.write('\t%s' % str(x))
			outfile.write('\n')

def write_stats(args, stats, species_info):
	# open output file
	outfile = open('%s/species_prevalence.txt' % args['outdir'], 'w')
	fields = ['mean_coverage', 'median_coverage', 'mean_abundance', 'median_abundance', 'prevalence']
	header = ['species_id', 'species_name']+fields
	outfile.write('\t'.join(header)+'\n')
	# order species ids by decreasing abundance
	sorted_species = [[x, y['prevalence']] for x, y in stats.items()]
	sorted_species.sort(key=lambda x: x[1], reverse=True)
	# write stats
	for species_id, prevalence in sorted_species:
		outfile.write(species_id)
		outfile.write('\t%s' % species_info[species_id])
		for field in fields:
			if field == 'prevalence':
				outfile.write('\t%s' % str(stats[species_id][field]))
			else:
				outfile.write('\t%s' % str(round(stats[species_id][field], 2)))
		outfile.write('\n')

def identify_samples(args):
	samples = []
	for sample in merge.load_samples(args):
		if not sample.paths['species']:
			sys.stderr.write("Warning: no species profile for %s\n" % sample.dir)
		elif sample.id in [s.id for s in samples]:
			sys.stderr.write("Warning: sample_id '%s' specified more than one time.\nSkipping: %s\n" % (sample.id, sample.dir))
		else:
			samples.append(sample)
	if len(samples)==0:
		sys.exit("\nError: no samples with species profiles\n")
	# select a subset of species to analyze
	if args['max_samples'] is not None and len(samples) > args['max_samples']:
		samples = samples[0:args['max_samples']]
	return samples

def run_pipeline(args):
	# list samples and species
	samples = identify_samples(args)
	species_info = species.read_annotations(args)
	# read in data & compute stats
	data = store_data(args, samples, species_info)
	stats = compute_stats(args, data)
	# write results
	write_abundance(args, samples, species_info, data)
	write_stats(args, stats, species_info)
