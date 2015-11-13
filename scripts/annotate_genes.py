#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.2'

import argparse, sys, os, gzip

def parse_arguments():
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__),
		description="""
			Transform a pan-genome matrix, where rows are genes, to a function matrix, where rows are functions. 
			Types of functions include KEGG Pathways, FIGFam gene families, Gene Ontology terms, or reactions from the Enzyme Commision database. 
			The abundances of genes (presence/absence or copy number) are summed by function_id. 
			Genes that have no annotated function are dropped.""")
	parser.add_argument('-i', dest='cnv_matrix', type=str, required=True,
						help='Gene CNV matrix. Expected file name: {species_id}.presabs or {species_id}.copynum')
	parser.add_argument('-o', dest='function_matrix', type=str, required=True,
						help='Function CNV matrix')
	parser.add_argument('-f', dest='ontology', choices=['kegg', 'figfams', 'go', 'ec'], required=True,
						help='kegg=KEGG pathways, figfams=FIGfams, go=Gene Ontology, ec=Enzyme Commission')
	parser.add_argument('-v', '--verbose', action='store_true', default=False)
	args = vars(parser.parse_args())
	check_args(args)
	args['cluster_id'] = os.path.basename(args['cnv_matrix']).split('.')[0]
	args['ref_db'] = '%s/ref_db/genome_clusters' % os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	return args

def check_args(args):
	if not os.path.isfile(args['cnv_matrix']):
		sys.exit("Input file does not exist!")
	if not os.path.isdir(os.path.dirname(os.path.abspath(args['function_matrix']))):
		sys.exit("Output directory does not exist!")

def read_gene_map(args):
	gene_to_functions = {}
	inpath = '%s/%s/functions.txt.gz' % (args['ref_db'], args['cluster_id'])
	infile = gzip.open(inpath)
	for index, line in enumerate(infile):
		gene_id, function_id, ontology = line.rstrip().split()
		if ontology in args['ontology']:
			if gene_id not in gene_to_functions:
				gene_to_functions[gene_id] = []
			gene_to_functions[gene_id].append(function_id)
	return gene_to_functions

def compute_abundances(args, gene_to_function):
	# open infile
	infile = open(args['cnv_matrix'])
	fields = next(infile).rstrip().split()
	# setup dict
	samples = fields[1:]
	functions = set([x for y in gene_to_function.values() for x in y])
	abundances = {}
	for sample in samples:
		abundances[sample] = {}
		for function_id in functions:
			abundances[sample][function_id] = 0.0
	# aggregate abundances
	for line in infile:
		record = line.rstrip().split()
		gene_id = record[0]
		if gene_id in gene_to_function:
			values = record[1:]
			for sample, value in zip(samples, values):
				for function_id in gene_to_function[gene_id]:
					abundances[sample][function_id] += float(value)
	return abundances

def write_results(args, gene_to_function, abundances):
	samples = abundances.keys()
	functions = set([x for y in gene_to_function.values() for x in y])
	outfile = open(args['function_matrix'], 'w')
	outfile.write('\t'.join(['function_id']+samples)+'\n')
	for function_id in functions:
		outfile.write(function_id)
		for sample_id in samples:
			outfile.write('\t%s' % abundances[sample_id][function_id])
		outfile.write('\n')

if __name__ == '__main__':
	
	args = parse_arguments()
	gene_to_function = read_gene_map(args)
	abundances = compute_abundances(args, gene_to_function)
	write_results(args, gene_to_function, abundances)



	




