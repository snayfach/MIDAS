#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip
from collections import defaultdict
from midas import utility
from midas.merge import merge

def read_gene_map(species_id, args):
	""" Map 99% centroids to gene_ids at lower level """
	gene_to_family = {}
	inpath = '%s/genome_clusters/%s/gene_family_map.txt.gz' % (args['db'], species_id)
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

def build_gene_matrices(species_id, samples, args):
	""" Compute gene copy numbers for samples """
	gene_to_family = read_gene_map(species_id, args)
	count_genes = len(gene_to_family.keys())
	count_genomes = len(set(['.'.join(x.split('.')[0:2]) for x in gene_to_family]))
	count_families = len(set(gene_to_family.values()))
	print("    %s genes from %s genomes" % (count_genes, count_genomes))
	print("    clustered into %s families at %s percent id" % (count_families, args['cluster_pid']))
	for sample in samples:
		sample.genes = {}
		for type in ['presabs', 'copynum', 'depth']:
			sample.genes[type] = defaultdict(float)
		inpath = '%s/genes/output/%s.genes.gz' % (sample.dir, species_id)
		for r in utility.parse_file(inpath):
			gene_id = gene_to_family[r['gene_id']]
			sample.genes['copynum'][gene_id] += float(r['normalized_coverage'])
			sample.genes['depth'][gene_id] += float(r['raw_coverage'])
	for sample in samples:
		for gene_id, copynum in sample.genes['copynum'].items():
			if copynum >= args['min_copy']: sample.genes['presabs'][gene_id] = 1
			else: sample.genes['presabs'][gene_id] = 0

def write_gene_matrices(species_id, samples, args):
	""" Compute pangenome matrices to file """
	# open outfiles
	outfiles = {}
	for type in ['presabs', 'copynum', 'depth']:
		outfiles[type] = open('%s/%s/%s.genes.%s' % (args['outdir'], species_id, species_id, type), 'w')
		outfiles[type].write('\t'.join(['gene_id'] + [s.id for s in samples])+'\n')
	# write values
	genes = sorted(samples[0].genes['depth'])
	for gene_id in genes:
		for type in ['presabs', 'copynum', 'depth']:
			outfiles[type].write(gene_id)
			for sample in samples:
				outfiles[type].write('\t%s' % str(sample.genes[type][gene_id]))
			outfiles[type].write('\n')
	for outfile in outfiles.values():
		outfile.close()

def write_gene_info(species_id, args):
	""" Write gene info file """
	gene_to_family = read_gene_map(species_id, args)
	outfile = open('%s/%s/%s.genes.info' % (args['outdir'], species_id, species_id), 'w')
	fields = ['gene_id', 'family_id', 'function_id', 'function_db']
	outfile.write('\t'.join(fields)+'\n')
	for ontology in ['figfam', 'kegg', 'go', 'ec']:
		gene_to_functions = read_function_map(args['db'], species_id, ontology)
		for gene_id, function_ids in gene_to_functions.items():
			family_id = gene_to_family[gene_id]
			for funtion_id in function_ids:
				outfile.write('\t'.join([gene_id, family_id, funtion_id, ontology])+'\n')

def run_pipeline(args):

	print("Identifying species")
	species = merge.select_species(args, type='genes')

	for sp in species:

		print "Merging: %s (id:%s) for %s samples" % (sp.consensus_name, sp.id, len(sp.samples))
		outdir = os.path.join(args['outdir'], sp.id)
		if not os.path.isdir(outdir): os.mkdir(outdir)
			
		print("  building pangenome matrices")
		build_gene_matrices(sp.id, sp.samples, args)
		write_gene_matrices(sp.id, sp.samples, args)
		
		print("  writing gene info file")
		write_gene_info(sp.id, args)

		print("  writing summary statistics")
		merge.write_summary_stats(sp.id, sp.samples, args, 'genes')

		print("")


