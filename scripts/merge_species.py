#!/usr/bin/env python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, sys
import argparse
import numpy as np
from phylo_cnv import species
from phylo_cnv import merge
from phylo_cnv import utility

def parse_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: merge_species.py [options]

Description: Merge species abundance files across samples
Input: list of sample directories
Output: relative abundance matrix, genome-coverage matrix, read-count matrix, species prevalence
""",
		epilog="""Examples:
1) provide list of paths to sample directories:
merge_species.py species -i /path/to/samples/sample_1,/path/to/samples/sample_2 -t list -o outbase

2) provide directory containing all samples:
merge_species.py species -i /path/to/samples -t dir -o outbase

3) provide file containing paths to sample directoriess:
merge_species.py species -i /path/to/samples/sample_paths.txt -t file -o outbase
""")
	parser.add_argument('-i', type=str, dest='input', required=True,
		help="""input to sample directories output by run_phylo_cnv.py species
can be a list of directories, a directory containing all samples, or a file with paths
see '-t' for details""")
	parser.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True,
		help="""'list': -i incdicates a comma-separated list of paths to sample directories
        example: /path/to/samples/sample_1,/path/to/samples/sample_2
'dir': -i incdicates a  directory containing all samples
       example: /path/to/samples
'file': -i incdicates a file containing paths to sample directories
	   example: /path/to/sample_paths.txt
""")
	parser.add_argument('-o', dest='outbase', type=str, required=True,
		help="basename for output files")
	parser.add_argument('-m', dest='min_cov', type=float, required=False, default=1.0,
		help="""minimum genome-coverage for estimating species prevalence (1.0)""")
	args = vars(parser.parse_args())
	args = utility.add_ref_db(args)
	check_args(args)
	return args

def print_arguments(args):
	print ("===========Parameters===========")
	print ("Script: merge_species.py")
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output basename: %s" % args['outbase'])
	print ("Minimum coverage for estimating prevalence: %s" % args['min_cov'])
	print ("")

def check_args(args):
	base_dir = os.path.dirname(os.path.abspath(args['outbase']))
	if not os.path.isdir(base_dir):
		sys.exit("Output directory of base name (-o) not found:\n%s" % base_dir)

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
		outfile = open('%s.%s' % (args['outbase'], field), 'w')
		outfile.write('\t'.join(['species_id']+[s.id for s in samples])+'\n')
		for species_id in species_info:
			outfile.write(species_id)
			for x in data[species_id][field]:
				outfile.write('\t%s' % str(x))
			outfile.write('\n')

def write_stats(args, stats, species_info):
	# open output file
	outfile = open('%s.species_prevalence' % args['outbase'], 'w')
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
	return samples

if __name__ == "__main__":

	args = parse_arguments()
	utility.print_copyright()
	print_arguments(args)
	
	# list samples and species
	samples = identify_samples(args)
	species_info = species.read_annotations(args)
	
	# read in data & compute stats
	data = store_data(args, samples, species_info)
	stats = compute_stats(args, data)

	# write results
	write_abundance(args, samples, species_info, data)
	write_stats(args, stats, species_info)

