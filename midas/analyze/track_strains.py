#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, itertools
from operator import itemgetter
from midas import utility
from midas.analyze import parse_snps

def id_markers(args):
	""" Pipeline for identifying marker alleles """

	# initialize input data
	species = parse_snps.Species(args['indir'])
	samples = parse_snps.fetch_samples(species, keep_samples=args['samples'])

	# open output file & write header
	outfile = open(args['out'], 'w')
	header = ['site_id', 'allele', 'count_samples'] + ['count_'+b for b in list('ATCG')]
	outfile.write('\t'.join(header)+'\n')
	
	# loop over sites & identify markers
	count_markers = 0
	sites = parse_snps.fetch_sites(species, samples)
	for index, site in enumerate(sites):
		
		# stop early
		if index >= args['max_sites']: break

		# count the occurence of 4 nucleotides at site across samples
		counts, total = count_alleles(site, args['samples'], args['min_freq'], args['min_reads'])

		# skip sites with non-discriminative allele(s)
		alleles = [(base, count) for base, count in counts.items() if count > 0]
		alleles.sort(key=itemgetter(1))
		if len(alleles) != 2: continue
		if alleles[0][1] > args['allele_prev']: continue

		# write identified discriminative allele
		count_markers += 1
		record = [site.id, alleles[0][0], total] + [counts[b] for b in list('ATCG')]
		outfile.write('\t'.join([str(_) for _ in record])+'\n')

	print("\n%s total disriminative alleles found" % count_markers)

def count_alleles(site, samples, min_freq, min_reads):
	""" compute the prevelance of the 4 nucleotides across samples """
	total = set([])
	groups = {'A':set([]), 'T':set([]), 'C':set([]), 'G':set([])}
	for sample in site.samples.values():
		if sample.depth == 0:
			continue
		if sample.freq >= min_freq and round(sample.freq * sample.depth) >= min_reads:
			groups[site.minor_allele].add(sample.id)
		if (1-sample.freq) >= min_freq and round((1-sample.freq) * sample.depth) >= min_reads:
			groups[site.major_allele].add(sample.id)
		total.add(sample.id)
	counts = dict([(allele, len(group)) for allele, group in groups.items()])
	return counts, len(total)

def track_markers(args):
	# initialize input data
	species = parse_snps.Species(args['indir'])
	samples = parse_snps.fetch_samples(species)
	species.paths['markers'] = args['markers']

	# open output file
	outfile = open(args['out'], 'w')
	header = ['sample1', 'sample2', 'count1', 'count2', 'count_both', 'count_either']
	outfile.write('\t'.join(header)+'\n')

	# determine marker alleles present in each sample
	print("Determining marker alleles present in each sample")
	call_markers(species, samples, args)

	# quantify marker allele sharing between samples
	print("Quantifying sharing of marker alleles between samples")
	allele_sharing(samples, outfile)

def call_markers(species, samples, args):
	""" determine if marker present in each sample """
	
	# open marker list
	markers = utility.parse_file(species.paths['markers'])
	marker = fetch_marker(markers) # dictionary for 1st marker allele
	if marker is None:
		sys.exit("\nError: no marker alleles found in file: %s\n" % species.paths['markers'])
	
	# init markers per sample
	for sample in samples.values():
		sample.markers = set([])

	# loop over sites
	sites = parse_snps.fetch_sites(species, samples)
	for index, site in enumerate(sites):
		
		# stop early
		if index >= args['max_sites']: break
		
		# skip sites not in marker list
		if (site.id != marker['site_id']):
			continue
			
		# determine if marker present in each sample
		for sample in site.samples.values():
		
			# skip samples without marker
			if sample.depth == 0:
				continue
			elif marker['allele'] == site.major_allele:
				sample.marker_freq = 1-sample.freq
			elif marker['allele'] == site.minor_allele:
				sample.marker_freq = sample.freq
			else:
				continue

			# record marker allele for sample_id
			sample.marker_count = round(sample.marker_freq * sample.depth)
			if (sample.marker_freq >= args['min_freq']
					and sample.marker_count >= args['min_reads']):
				sample.markers.add(site.id)
				
		# fetch next marker allele
		marker = fetch_marker(markers)
		if marker is None: break

def fetch_marker(markers):
	""" Fetch next marker allele from file """
	try:
		marker = next(markers)
		return marker
	except StopIteration:
		return None

def allele_sharing(samples, outfile):
	""" Compute sharing between sample pairs """
	for index, pair in enumerate(itertools.combinations(samples, r=2)):
		if not index % 500: print("%s sample pairs processed" % index)
		sample1, sample2 = pair
		alleles1 = samples[sample1].markers
		alleles2 = samples[sample2].markers
		count1 = len(alleles1)
		count2 = len(alleles2)
		count_both = len(set(alleles1 & alleles2))
		count_either = len(set(alleles1 | alleles2))
		record = [sample1, sample2, count1, count2, count_both, count_either]
		outfile.write('\t'.join([str(_) for _ in record])+'\n')


