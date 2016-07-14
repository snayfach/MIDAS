#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os
from operator import itemgetter
from midas import parse

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: 
  1) compute the prevalence of the 4 nucleotides for a species across a set of specified sample ids
  2) select a subset of discriminative alleles at bi-allelic sites
  3) write these alleles to specified output file
  
Output fields:
  site_id - site identifier (scaffold_id|position)
  allele - nucleotide (A, T, C, or G)
  count_samples - number of samples with non-zero coverage at position
  count_A - number of samples with A
  count_T - number of samples with T
  count_C - number of samples with C
  count_G - number of samples with G

Usage: id_marker_alleles.py indir outpath sample_to_group [options]
""",
	epilog="""""")
	parser.add_argument('indir', type=str,
		help="""path to input directory
contains files: [snps_alt_allele.txt, snps_info.txt, snps_summary.txt
                 snps_depth.txt, snps_ref_freq.txt]""")
	parser.add_argument('outpath', type=str,
		help="""path to output file""")
	parser.add_argument('sample_to_group', type=str,
		help="""tab delimited file with no header and two columns: sample_id and group_id
sample_id denotes sample ids present in INDIR
group_id denotes non-independent sample_ids (e.g. technical, biological replicates)
marker alleles are chosen such that they occur once per independent sample_id (see below)
example sample_to_group file:
sample1 subject1
sample2 subject1
sample3 subject2
sample4 subject2
""")
	presabs = parser.add_argument_group("Options for determining allele presence/absence")
	presabs.add_argument('--min_freq', type=float, metavar='FLOAT', default=0.10,
		help="""minimum frequency (proportion of reads) for determining the presence/absence
of a nucleotide at a genomic site (0.10)""")
	presabs.add_argument('--min_reads', type=int, metavar='INT', default=3,
		help="""minimum number of reads for determining the presence/absence
of a nucleotide at a genomic site (3)""")

	sites = parser.add_argument_group("Options for selecting discriminative sites")
	sites.add_argument('--max_groups', type=int, metavar='INT', default=1,
		help="""discard alleles found in > MAX_GROUPS (1)
a value of one indicates that marker alleles only occur once per independent sample""")
	sites.add_argument('--max_sites', type=int, metavar='INT',
		help="""maximum number of genomic sites to process (use all)
useful for quick tests""")
	args = vars(parser.parse_args())
	args['samples'] = dict([_.rstrip().split() for _ in open(args['sample_to_group'])])
	return args

def allele_counts(site, samples, min_freq, min_reads):
	""" compute the prevelance of the 4 nucleotides for site across samples """
	total = set([]) # number of samples with non-zero coverage for site
	groups = {'A':set([]), 'T':set([]), 'C':set([]), 'G':set([])}
	for sample_id, group_id in samples.items(): # fetch site info for each sample
		depth = int(site.depth[sample_id])
		if depth == 0: # skip samples with no reads
			continue
		ref_freq = float(site.ref_freq[sample_id]) # count reference allele
		if ref_freq >= min_freq and round(ref_freq * depth) >= min_reads:
			ref_allele = site.info['ref_allele']
			groups[ref_allele].add(group_id)
		alt_freq = 1-ref_freq # count alternate allele
		if alt_freq >= min_freq and round(alt_freq * depth) >= min_reads:
			alt_allele = site.alt_allele[sample_id]
			groups[alt_allele].add(group_id)
		total.add(group_id) # keep track of number of samples with data
	counts = dict([(allele, len(group)) for allele, group in groups.items()])
	return counts, len(total)

def init_paths(args):
	""" fetch paths to input files """
	paths = {}
	exts = ['alt_allele', 'info', 'summary', 'depth', 'ref_freq']
	for ext in exts:
		inpath = '%s/snps_%s.txt' % (args['indir'], ext)
		if os.path.isfile(inpath):
			paths[ext] = inpath
		else:
			sys.exit("Input file does not exist: %s" % inpath)
	return paths

def check_sample_ids(args, paths):
	""" make sure that specified sample ids are present in input """
	samples = [_.split()[0] for _ in open(paths['summary'])][1:]
	for sample_id in args['samples'].copy():
		if sample_id not in samples:
			del args['samples'][sample_id]
	if len(args['samples']) == 0:
		sys.exit("No samples remain")

if __name__ == "__main__":

	args = parse_arguments() # dictionary of args
	
	paths = init_paths(args) # dictionary of paths
	
	check_sample_ids(args, paths) # make sure that specified sample ids are present in input
	
	# open output file & write header
	bases = ['A', 'T', 'C', 'G']
	outfile = open(args['outpath'], 'w')
	header = ['site_id', 'allele', 'count_samples'] + ['count_'+b for b in bases]
	outfile.write('\t'.join(header)+'\n')

	count_alleles = 0
	for index, site in enumerate(parse.parse_sites(args, paths)):
		
		# record progress
		if not index % 100000: print("%s sites processed" % index)
		
		# count the occurence of 4 nucleotides at site across samples
		counts, total = allele_counts(site, args['samples'], args['min_freq'], args['min_reads'])
		
		# skip sites with non-discriminative allele(s)
		alleles = [(base, count) for base, count in counts.items() if count > 0]
		alleles.sort(key=itemgetter(1))
		if len(alleles) != 2:
			continue
		elif alleles[0][1] > args['max_groups']:
			continue

		# write identified discriminative allele
		count_alleles += 1
		record = [site.id, alleles[0][0], total] + [counts[b] for b in bases]
		outfile.write('\t'.join([str(_) for _ in record])+'\n')

	print("%s total disriminative alleles found" % count_alleles)





