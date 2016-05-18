#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip
from collections import defaultdict
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import merge, utility

def identify_samples(args):
	""" Identify samples with sufficient depth for target genome-cluster """
	samples = []
	for sample in merge.load_samples(args):
		if not sample.paths['genes']:
			sys.stderr.write("Warning: no genes output for sample: %s\n" % sample.dir)
			continue
		stats = merge.read_stats(sample.paths['genes'])
		if args['species_id'] not in stats:
			continue
		elif float(stats[args['species_id']]['phyeco_coverage']) < args['marker_coverage']:
			continue
		elif float(stats[args['species_id']]['mean_coverage']) < args['gene_coverage']:
			continue
		else:
			samples.append(sample)
		if args['max_samples'] and len(samples) >= args['max_samples']:
			break
	if len(samples) == 0:
		sys.exit("\nError: no samples with sufficient coverage for species_id.\nTry running with more lenient parameters")
	else:
		print("  %s samples with species" % len(samples))
		return samples

def read_gene_map(args):
	""" Map 99% centroids to gene_ids at lower level """
	gene_to_family = {}
	inpath = '%s/genome_clusters/%s/gene_family_map.txt.gz' % (args['db'], args['species_id'])
	infile = gzip.open(inpath)
	fields = next(infile).rstrip().split()
	for line in infile:
		values = line.rstrip().split()
		map = dict([(f,v) for f,v in zip(fields, values)])
		gene_to_family[map['99']] = map[args['cluster_pid']]
	return gene_to_family

def read_function_map(ref_db, species_id, ontology):
	""" Map gene ids to functions for given ontology """
	gene_to_functions = {}
	inpath = '%s/genome_clusters/%s/functions.txt.gz' % (ref_db, species_id)
	infile = gzip.open(inpath)
	for index, line in enumerate(infile):
		gene_id, function_id, ont = line.rstrip().split()
		if ont == ontology:
			if gene_id not in gene_to_functions:
				gene_to_functions[gene_id] = []
			gene_to_functions[gene_id].append(function_id)
	return gene_to_functions

def build_gene_matrices(samples, args):
	""" Compute gene copy numbers for samples """
	gene_to_family = read_gene_map(args)
	count_genes = len(gene_to_family.keys())
	count_genomes = len(set(['.'.join(x.split('.')[0:2]) for x in gene_to_family]))
	count_families = len(set(gene_to_family.values()))
	print("  %s genes from %s genomes" % (count_genes, count_genomes))
	print("  clustered into %s families at %s percent id" % (count_families, args['cluster_pid']))
	for sample in samples:
		sample.genes = {}
		for type in ['presabs', 'copynum', 'depth']:
			sample.genes[type] = defaultdict(float)
		inpath = '%s/genes/coverage/%s.cov.gz' % (sample.dir, args['species_id'])
		for r in utility.parse_file(inpath):
			gene_id = gene_to_family[r['gene_id']]
			sample.genes['copynum'][gene_id] += float(r['normalized_coverage'])
			sample.genes['depth'][gene_id] += float(r['raw_coverage'])
	for sample in samples:
		for gene_id, copynum in sample.genes['copynum'].items():
			if copynum >= args['min_copy']: sample.genes['presabs'][gene_id] = 1
			else: sample.genes['presabs'][gene_id] = 0

def write_gene_matrices(samples, args):
	""" Compute pangenome matrices to file """
	# open outfiles
	outfiles = {}
	for type in ['presabs', 'copynum', 'depth']:
		outfiles[type] = open('%s/%s.genes.pangenome.%s' % (args['outdir'], args['species_id'], type), 'w')
		outfiles[type].write('\t'.join(['gene_id'] + [s.id for s in samples])+'\n')
	# write values
	genes = sorted(samples[0].genes['depth'])
	for gene_id in genes:
		for type in ['presabs', 'copynum', 'depth']:
			outfiles[type].write(gene_id)
			for sample in samples:
				outfiles[type].write('\t%s' % str(sample.genes[type][gene_id]))
			outfiles[type].write('\n')

def build_function_matrices(samples, ontology, args):
	""" Compute function copy numbers for samples """
	gene_to_functions = read_function_map(args['db'], args['species_id'], ontology)
	count_genes = len(gene_to_functions.keys())
	count_functions = len(set([i for s in gene_to_functions.values() for i in s]))
	print("  %s genes mapped to %s functions from %s ontology" % (count_genes, count_functions, ontology))
	for sample in samples:
		sample.functions = {}
		for type in ['presabs', 'copynum', 'depth']:
			sample.functions[type] = defaultdict(float)
			for gene_id, value in sample.genes[type].items():
				if gene_id in gene_to_functions:
					for function_id in gene_to_functions[gene_id]:
						sample.functions[type][function_id] += value

def write_function_matrices(samples, ontology, args):
	""" Compute pangenome matrices to file """
	# open outfiles
	outfiles = {}
	for type in ['presabs', 'copynum', 'depth']:
		outfiles[type] = open('%s/%s.genes.%s.%s' % (args['outdir'], args['species_id'], ontology, type), 'w')
		outfiles[type].write('\t'.join(['gene_id'] + [s.id for s in samples])+'\n')
	# write values
	functions = sorted(samples[0].functions['depth'])
	for function_id in functions:
		for type in ['presabs', 'copynum', 'depth']:
			outfiles[type].write(function_id)
			for sample in samples:
				outfiles[type].write('\t%s' % str(sample.functions[type][function_id]))
			outfiles[type].write('\n')

def write_summary_stats(samples, args):
	outfile = open('%s/%s.genes.summary' % (args['outdir'], args['species_id']), 'w')
	fields = ['pangenome_size', 'covered_genes', 'mean_coverage', 'phyeco_coverage']
	outfile.write('\t'.join(['species_id', 'sample_id']+fields)+'\n')
	for sample in samples:
		stats = merge.read_stats(sample.paths['genes'])
		outfile.write('\t'.join([args['species_id'], sample.id]))
		for field in fields:
			outfile.write('\t%s' % str(stats[args['species_id']][field]))
		outfile.write('\n')

def run_pipeline(args):

	print("Identifying samples with species")
	samples = identify_samples(args)
		
	print("Building pangenome matrices")
	build_gene_matrices(samples, args)
	write_gene_matrices(samples, args)
	
	print("Building function matrices")
	for ontology in ['figfam', 'kegg', 'go', 'ec']:
		if args['no_functions']: continue
		build_function_matrices(samples, ontology, args)
		write_function_matrices(samples, ontology, args)

	print("Writing summary statistics")
	write_summary_stats(samples, args)


