#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip
from collections import defaultdict
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import merge, utility

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: merge_genes.py [options]

Description: merge results from pan-genome profiling across samples
Input: list of sample directories
Output: pan-genome copy-number matrix, presence/absence matrix, and read-depth matrix
        matrixes also created for KEGG, FIGfams, Gene Ontology, and Enzyme Comission (E.C.)
""",
		epilog="""Examples:
1) Merge results for species 57955. Provide list of paths to sample directories:
merge_genes.py -s 57955 -o outdir/57955 -i sample_1,sample_2 -t list

2) Build matrix for pan-genome genes at lower percent id threshold:
merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --cluster_pid 85

3) Exclude low-coverage samples in output matrix:
merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --marker_coverage 5.0

4) Use lenient threshold for determining gene presence-absence:
merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --min_copy 0.1

5) Just write pan-genome matrices; do not write results for KEGG, FIGfams, GO, or EC:
merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --no_functions
""")
	
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""input to sample directories output by run_phylo_cnv.py genes
see '-t' for details""")
	io.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True,
		help="""'list': -i is a comma-separated list of paths to sample directories (ex: /sample1,/sample2)
'dir': -i is a  directory containing all samples (ex: /samples_dir)
'file': -i is a file containing paths to sample directories (ex: sample_paths.txt)
""")
	io.add_argument('-s', dest='species_id', type=str, required=True,
		help="""species identifier
a list of prevalent species can be obtained by running 'merge_species.py'
a map of species ids to species names can be found in 'ref_db/annotations.txt'""")
	io.add_argument('-o', type=str, dest='outdir', required=True,
		help="""output directory""")
	io.add_argument('--no_functions', action='store_true', default=False,
		help="""do not write function matrices to OUTDIR""")

	sample = parser.add_argument_group('Sample filters (select subset of samples from INPUT)')
	sample.add_argument('--marker_coverage', type=float, default=1.0, metavar='FLOAT',
		help="""minimum coverage per sample across 15 phylogenetic marker genes (1.0)""")
	sample.add_argument('--gene_coverage', type=float, default=1.0, metavar='FLOAT',
		help="""minimum coverage per sample across all genes with non-zero coverage (1.0)""")
	sample.add_argument('--max_samples', type=int, metavar='INT',
		help="""maximum number of samples to process. useful for testing (use all)""")

	gene = parser.add_argument_group('Presence/Absence')
	gene.add_argument('--min_copy', type=float, default=0.35, metavar='FLOAT',
		help="""genes >= MIN_COPY are classified as present
genes < MIN_COPY are classified as absent (0.35)""")
	gene.add_argument('--cluster_pid', type=str, dest='cluster_pid', default='95', choices=['75', '80', '85', '90', '95', '99'],
		help="""gene family percent identity
small values: fewer, larger gene families
large values: more, smaller gene families (95)""")
	
	args = vars(parser.parse_args())
	args = utility.add_ref_db(args)
	
	return args

def print_arguments(args):
	print ("===========Parameters===========")
	print ("Script: merge_genes.py")
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Species identifier: %s" % args['species_id'])
	print ("Sample selection criteria:")
	if args['marker_coverage']:
		print ("  >=%s average coverage across 15 universal-single-copy genes" % args['marker_coverage'])
	if args['gene_coverage']:
		print ("  >=%s average coverage across all genes with non-zero coverage" % args['gene_coverage'])
	if args['max_samples']:
		print ("  analyze up to %s samples" % args['max_samples'])
	print ("Gene quantification criterea:")
	print ("  present (1): genes with copy number >=%s" % args['min_copy'])
	print ("  absent (0): genes with copy number <%s" % args['min_copy'])
	print ("  cluster genes at %s percent identity" % args['cluster_pid'])
	print ("")

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
	outfile = open('%s/%s.genes.summary_stats' % (args['outdir'], args['species_id']), 'w')
	fields = ['pangenome_size', 'covered_genes', 'phyeco_coverage', 'mean_coverage']
	outfile.write('\t'.join(['species_id', 'sample_id']+fields)+'\n')
	for sample in samples:
		stats = merge.read_stats(sample.paths['genes'])
		outfile.write('\t'.join([args['species_id'], sample.id]))
		for field in fields:
			outfile.write('\t%s' % str(stats[args['species_id']][field]))
		outfile.write('\n')

if __name__ == '__main__':

	args = parse_arguments()
	utility.print_copyright()
	print_arguments(args)
	if not os.path.isdir(args['outdir']): os.mkdir(args['outdir'])

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


