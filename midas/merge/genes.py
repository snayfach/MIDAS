#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip
from collections import defaultdict
from midas import utility
from midas.merge import merge

def build_gene_matrices(sp, min_copy):
	""" Compute gene copy numbers for samples """
	for sample in sp.samples:
		sample.genes = {}
		for field, dtype in [('presabs',float), ('copynum',float), ('depth',float), ('reads',int)]:
			sample.genes[field] = defaultdict(dtype)
		inpath = '%s/genes/output/%s.genes.gz' % (sample.dir, sp.id)
		for r in utility.parse_file(inpath):
			if 'ref_id' in r: r['gene_id'] = r['ref_id'] # fix old fields if present
			if 'normalized_coverage' in r: r['copy_number'] = r['normalized_coverage'] 
			if 'raw_coverage' in r: r['coverage'] = r['raw_coverage']
			gene_id = sp.map[r['gene_id']]
			sample.genes['copynum'][gene_id] += float(r['copy_number'])
			sample.genes['depth'][gene_id] += float(r['coverage'])
			sample.genes['reads'][gene_id] += int(r['count_reads']) if 'count_reads' in r else 0
	for sample in sp.samples:
		for gene_id, copynum in sample.genes['copynum'].items():
			if copynum >= min_copy: sample.genes['presabs'][gene_id] = 1
			else: sample.genes['presabs'][gene_id] = 0

def write_gene_matrices(sp):
	""" Compute pangenome matrices to file """
	# open outfiles
	outfiles = {}
	for type in ['presabs', 'copynum', 'depth', 'reads']:
		outfiles[type] = open('%s/genes_%s.txt' % (sp.dir, type), 'w')
		outfiles[type].write('\t'.join(['gene_id'] + [s.id for s in sp.samples])+'\n')
	# write values
	genes = sorted(sp.samples[0].genes['depth'])
	for gene_id in genes:
		for type in ['presabs', 'copynum', 'depth', 'reads']:
			outfiles[type].write(gene_id)
			for sample in sp.samples:
				outfiles[type].write('\t%s' % str(sample.genes[type][gene_id]))
			outfiles[type].write('\n')
	for outfile in outfiles.values():
		outfile.close()

def write_readme(args, sp):
	outfile = open('%s/readme.txt' % sp.dir, 'w')
	outfile.write("""
Description of output files and file formats from 'merge_midas.py genes'

Output files
############
genes_depth.txt  
  average-read depth of each gene per sample
genes_copynum.txt
  copy-number of each gene per sample
  estimated by dividing the read-depth of a gene by the median read-depth of 15 universal single copy genes
genes_presabs.txt  
  the presence (1) or absence (0) of each gene per sample
  estimated by applying a threshold to gene copy-number values
genes_reads.txt
  number of reads mapped to each gene per sample
genes_summary.txt
  alignment summary statistics per sample

Output formats
############
genes_depth.txt, genes_copynum.txt, genes_presabs.txt, genes_reads.txt
  tab-delimited matrix files
  field names are sample ids
  row names are gene ids
genes_summary.txt
  sample_id: sample identifier
  pangenome_size: number of non-redundant genes in reference pan-genome
  covered_genes: number of genes with at least 1 mapped read
  fraction_covered: proportion of genes with at least 1 mapped read
  mean_coverage: average read-depth across genes with at least 1 mapped read
  marker_coverage: median read-depth across 15 universal single copy genes
  aligned_reads: number of reads that aligned to pangenome
  mapped_reads: number of aligned reads after applying filters for mapping quality, base quality, alignment fraction, and percent identity

Additional information for species can be found in the reference database:
 %s/pan_genomes/%s
""" % (args['db'], sp.id) )
	outfile.close()

def read_cluster_map(sp, db, pid):
	sp.map = {}
	for ext in ['', '.gz']:
		path = '/'.join([db, 'pan_genomes', sp.id, 'gene_info.txt%s' % ext])
		if os.path.isfile(path):
			sp.gene_info = path
	for r in utility.parse_file(sp.gene_info):
		sp.map[r['centroid_99']] =  r['centroid_%s' % pid]

def run_pipeline(args):

	print("Identifying species and samples")
	species_list = merge.select_species(args, dtype='genes')
	for species in species_list:
		print("  %s" % species.id)
		print("    count genomes: %s" % species.info['count_genomes'])
		print("    count samples: %s" % len(species.samples))
		
	print("\nMerging genes")
	for species in species_list:
			
		print("  %s" % species.id)
		species.dir = os.path.join(args['outdir'], species.id)
		if not os.path.isdir(species.dir): os.mkdir(species.dir)
		read_cluster_map(species, args['db'], args['cluster_pid'])
			
		print("    building pangenome matrices")
		build_gene_matrices(species, min_copy=args['min_copy'])
		write_gene_matrices(species)

		print("    writing summary statistics")
		species.write_sample_info(dtype='genes', outdir=args['outdir'])

		write_readme(args, species)
		
		print("    done!")


