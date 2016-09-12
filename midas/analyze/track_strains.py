
import argparse, sys, os, platform, itertools
from operator import itemgetter
from midas import utility
from midas.analyze import snp_matrix

def allele_counts(site, samples, min_freq, min_reads):
	""" compute the prevelance of the 4 nucleotides for site across samples """
	total = set([]) # number of samples with non-zero coverage for site
	groups = {'A':set([]), 'T':set([]), 'C':set([]), 'G':set([])}
	sample_site = site.sample_values()
	for sample_id, group_id in samples.items(): # fetch site info for each sample
		if sample_id not in sample_site: continue
		depth = int(sample_site[sample_id]['depth'])
		if depth == 0: # skip samples with no reads
			continue
		ref_freq = float(sample_site[sample_id]['ref_freq']) # count reference allele
		if ref_freq >= min_freq and round(ref_freq * depth) >= min_reads:
			groups[site.ref_allele].add(group_id)
		alt_freq = 1-ref_freq # count alternate allele
		if alt_freq >= min_freq and round(alt_freq * depth) >= min_reads:
			alt_allele = sample_site[sample_id]['alt_allele']
			groups[alt_allele].add(group_id)
		total.add(group_id) # keep track of number of samples with data
	counts = dict([(allele, len(group)) for allele, group in groups.items()])
	return counts, len(total)

def id_markers(args):
	# open output file & write header
	bases = ['A', 'T', 'C', 'G']
	outfile = open(args['out'], 'w')
	header = ['site_id', 'allele', 'count_samples'] + ['count_'+b for b in bases]
	outfile.write('\t'.join(header)+'\n')
	count_alleles = 0
	for index, site in enumerate(snp_matrix.parse_sites(args['indir'])):
		# record progress
		if not index % 100000: print("%s sites processed" % index)
		if args['max_sites'] and index >= args['max_sites']: break
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

def call_markers(args):
	""" determine if marker present in each sample """
	markers = utility.parse_file(args['markers']) # generator for marker alleles file
	marker = fetch_marker(markers) # dictionary for 1st marker allele
	samples = dict([(_, {}) for _ in snp_matrix.list_samples(args['indir'], args['max_samples'])]) # dic of marker alleles
	for index, site in enumerate(snp_matrix.parse_sites(args['indir'])):
		if not index % 100000: print("%s sites processed" % index)
		if args['max_sites'] is not None and index >= args['max_sites']: break
		# skip non-discriminative sites
		if (site.ref_id != marker['ref_id']
				or site.ref_pos < marker['ref_pos']):
			continue
		# determine if marker present in each sample
		sample_site = site.sample_values()
		for sample_id in samples:
			# skip samples without marker
			depth = int(sample_site[sample_id]['depth'])
			if depth == 0:
				continue
			elif marker['allele'] == site.ref_allele:
				freq = float(sample_site[sample_id]['ref_freq'])
			elif marker['allele'] == sample_site[sample_id]['alt_allele']:
				freq = 1-float(sample_site[sample_id]['ref_freq'])
			else:
				continue
			# record marker allele for sample_id
			if (freq >= args['min_freq']
					and round(freq * depth) >= args['min_reads']):
				samples[sample_id][site.id] = freq
		# fetch next marker allele
		marker = fetch_marker(markers)
		if marker is None: break # stop when there are no more markers
	return samples

def fetch_marker(markers):
	try:
		marker = next(markers)
		marker['ref_id'] = marker['site_id'].rsplit('|', 1)[0]
		marker['ref_pos'] = int(marker['site_id'].rsplit('|')[2])
		return marker
	except StopIteration:
		return None

def allele_sharing(x, y):
	a = len(x)
	b = len(y)
	i = len(set(x.keys()) & set(y.keys()))
	u = len(set(x.keys()) | set(y.keys()))
	return a, b, i, u

def track_markers(args):
	# open output file
	outfile = open(args['out'], 'w')
	header = ['sample1', 'sample2', 'count1', 'count2', 'count_both', 'count_either']
	outfile.write('\t'.join(header)+'\n')
	# determine marker alleles present in each sample
	print("Determining marker alleles present in each sample")
	samples = call_markers(args)
	# quantify marker allele sharing
	print("Quantifying sharing of marker alleles between samples")
	for index, pair in enumerate(itertools.combinations(samples, r=2)):
		if not index % 500: print("%s sample pairs processed" % index)
		sample1, sample2 = pair
		count1, count2, count_both, count_either = allele_sharing(samples[sample1], samples[sample2])
		record = [sample1, sample2, count1, count2, count_both, count_either]
		outfile.write('\t'.join([str(_) for _ in record])+'\n')

