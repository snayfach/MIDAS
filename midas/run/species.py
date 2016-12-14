#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, subprocess, Bio.SeqIO
from time import time
from midas import utility
from operator import itemgetter

def read_annotations(args):
	info = {}
	inpath = '%s/species_info.txt' % args['db']
	for r in utility.parse_file(inpath):
		info[r['species_id']] = r
	return info

def read_marker_info(args):
	""" Read info for marker genes from phyeco.fa """
	info = {}
	for seq in Bio.SeqIO.parse('%s/marker_genes/phyeco.fa' % args['db'], 'fasta'):
		info[seq.id] = None
	for r in utility.parse_file('%s/marker_genes/phyeco.map' % args['db']):
		if r['gene_id'] in info:
			info[r['gene_id']] = r
	return info

def map_reads_hsblast(args):
	""" Use hs-blastn to map reads in fasta file to marker database """
	# stream sequences
	command = 'python %s' % args['stream_seqs']
	command += ' -i %s' % args['m1'] # 1st mate
	if args['m2']: command += ',%s' % args['m2'] # 2nd mate
	if args['max_reads']: command += ' -n %s' % args['max_reads'] # number of reads
	if args['read_length']: command += ' -l %s' % args['read_length'] # read length
	command += ' 2> %s/species/temp/read_count.txt' % args['outdir'] # tmpfile to store # of reads, bp sampled
	# hs-blastn
	command += ' | %s align' % args['hs-blastn']
	command += ' -word_size %s' % args['word_size']
	command += ' -query /dev/stdin'
	command += ' -db %s/marker_genes/phyeco.fa' % args['db']
	command += ' -outfmt 6'
	command += ' -num_threads %s' % args['threads']
	command += ' -out %s/species/temp/alignments.m8' % args['outdir']
	command += ' -evalue 1e-3'
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)

def parse_blast(inpath):
	""" Yield formatted record from BLAST m8 file """
	formats = [str,str,float,int,float,float,float,float,float,float,float,float]
	fields = ['query','target','pid','aln','mis','gaps','qstart','qend','tstart','tend','evalue','score']
	for line in open(inpath):
		values = line.rstrip().split()
		yield dict([(field, format(value)) for field, format, value in zip(fields, formats, values)])

def query_coverage(aln):
	""" Compute alignment coverage of query """
	qlen = aln['query'].split('_')[-1] # get qlen from sequence header
	return float(aln['aln'])/int(qlen)

def find_best_hits(args, marker_info):
	""" Find top scoring alignment for each read """
	best_hits = {}
	marker_cutoffs = get_markers(args)
	i = 0
	qcovs = []
	for aln in parse_blast('%s/species/temp/alignments.m8' % args['outdir']):
		i += 1
		marker_id = marker_info[aln['target']]['marker_id'] # get gene family from marker_info
		cutoff = args['mapid'] if args['mapid'] else marker_cutoffs[marker_id]
		if aln['pid'] < cutoff: # does not meet marker cutoff
			continue
		elif query_coverage(aln) < args['aln_cov']: # filter local alignments
			continue
		elif aln['query'] not in best_hits: # record aln
			best_hits[aln['query']] = [aln]
		elif best_hits[aln['query']][0]['score'] == aln['score']: # add aln
			best_hits[aln['query']] += [aln]
		elif best_hits[aln['query']][0]['score'] < aln['score']: # update aln
			best_hits[aln['query']] = [aln]
	print("  total alignments: %s" % i)
	return best_hits.values()

def assign_unique(args, alns, species_info, marker_info):
	""" Count the number of uniquely mapped reads to each genome species """
	unique_alns = dict([(_,[]) for _ in species_info])
	unique = 0
	non_unique = 0
	for aln in alns:
		if len(aln) == 1:
			unique += 1
			#species_id = aln[0]['target'].split('_')[0]
			species_id = marker_info[aln[0]['target']]['species_id']
			unique_alns[species_id].append(aln[0])
		else:
			non_unique += 1
	print("  uniquely mapped reads: %s" % unique)
	print("  ambiguously mapped reads: %s" % non_unique)
	return unique_alns

def assign_non_unique(args, alns, unique_alns, marker_info):
	""" Probabalistically assign ambiguously mapped reads """
	import numpy as np
	import random
	total_alns = unique_alns.copy()
	for aln in alns:
		if len(aln) > 1:
			species_ids = [marker_info[_['target']]['species_id'] for _ in aln]
			counts = [len(unique_alns[_]) for _ in species_ids]
			if sum(counts) == 0:
				species_id = random.sample(species_ids, 1)[0]
			else:
				probs = [float(count)/sum(counts) for count in counts]
				species_id = np.random.choice(species_ids, 1, p=probs)[0]
			total_alns[species_id].append(aln[species_ids.index(species_id)])
	return total_alns

def get_markers(args):
	""" Read in optimal mapping parameters for marker genes; override if user has provided cutoff """
	marker_cutoffs = {}
	inpath = '/'.join([args['db'], 'marker_genes/phyeco.mapping_cutoffs'])
	if not os.path.isfile(inpath): sys.exit("File not found: %s" % inpath)
	for line in open(inpath):
		marker_id, min_pid = line.rstrip().split()
		if args['mapid']:
			marker_cutoffs[marker_id] = args['mapid']
		else:
			marker_cutoffs[marker_id] = float(min_pid)
	return marker_cutoffs

def read_gene_lengths(args, species_info, marker_info):
	""" Read in total gene length per species_id """
	total_gene_length = dict([(_,0) for _ in species_info])
	for r in marker_info.values():
		total_gene_length[r['species_id']] += int(r['gene_length'])
	return total_gene_length

def normalize_counts(species_alns, total_gene_length):
	""" Normalize counts by gene length and sum contrain """
	# norm by gene length, compute cov
	species_abundance = {}
	total_cov = 0.0
	for species_id, alns in species_alns.items():
		species_abundance[species_id] = {}
		# compute coverage
		if len(alns) > 0:
			bp = sum([aln['aln'] for aln in alns])
			cov = float(bp)/total_gene_length[species_id]
		else:
			cov = 0.0
		# store results
		species_abundance[species_id] = {'cov':cov, 'count':len(alns)}
		total_cov += cov
	# compute relative abundance
	total_cov = sum([_['cov'] for _ in species_abundance.values()])
	for species_id in species_abundance.keys():
		cov = species_abundance[species_id]['cov']
		species_abundance[species_id]['rel_abun'] = cov/total_cov if total_cov > 0 else 0
	print("  total marker-gene coverage: %s" % round(total_cov, 3))
	return species_abundance

def write_abundance(outdir, species_abundance, annotations):
	""" Write species results to specified output file """
	outpath = '%s/species/species_profile.txt' % outdir
	outfile = open(outpath, 'w')
	fields = ['species_id', 'count_reads', 'coverage', 'relative_abundance']
	outfile.write('\t'.join(fields)+'\n')
	species_ids =  sorted([(x, y['count']) for x, y in species_abundance.items()], key=itemgetter(1), reverse=True)
	for species_id, count_reads in species_ids:
		values = species_abundance[species_id]
		record = [species_id, values['count'], values['cov'], values['rel_abun']]
		outfile.write('\t'.join([str(x) for x in record])+'\n')

def read_abundance(inpath):
	""" Parse species abundance file """
	if not os.path.isfile(inpath):
		sys.exit("\nCould not locate species profile: %s\nTry rerunning with run_midas.py species" % inpath)
	abun = {}
	for rec in utility.parse_file(inpath):
		# format record
		if 'species_id' in rec: rec['species_id'] = rec['species_id']
		if 'count_reads' in rec: rec['count_reads'] = int(rec['count_reads'])
		if 'coverage' in rec: rec['coverage'] = float(rec['coverage'])
		if 'relative_abundance' in rec: rec['relative_abundance'] = float(rec['relative_abundance'])
		abun[rec['species_id']] = rec
	return abun

def select_species(args):
	""" Select genome species to map to """
	import operator
	species_sets = {}
	# read in species abundance if necessary
	if any([args['species_topn'], args['species_cov']]):
		species_abundance = read_abundance('%s/species/species_profile.txt' % args['outdir'])
		# user specifed a coverage threshold
		if args['species_cov']:
			species_sets['species_cov'] = set([])
			for species_id, values in species_abundance.items():
				if values['coverage'] >= args['species_cov']:
					species_sets['species_cov'].add(species_id)
		# user specifed topn genome-species
		if args['species_topn']:
			species_sets['species_topn'] = set([])
			species_abundance = [(i,d['relative_abundance']) for i,d in species_abundance.items()]
			sorted_abundance = sorted(species_abundance, key=operator.itemgetter(1), reverse=True)
			for species_id, rel_abun in sorted_abundance[0:args['species_topn']]:
				species_sets['species_topn'].add(species_id)
	# user specified a list of one or more genome-species
	if args['species_id']:
		species_sets['species_id'] = set([])
		for species_id in args['species_id']:
			species_sets['species_id'].add(species_id)
	# intersect sets of genome-species
	my_species = list(set.intersection(*species_sets.values()))
	# optionally remove bad species_ids
	inpath = '/'.join([args['db'], 'exclude.txt'])
	if os.path.isfile(inpath):
		for line in open(inpath):
			try: my_species.remove(line.rstrip())
			except: pass
	# check that at least one genome-species was selected
	if len(my_species) == 0:
		sys.exit("\nError: no species sastisfied your selection criteria. \n")
	return my_species
	
def scale_depth(species_abundance, mate1, mate2, file_type, max_reads):
	""" Scale read depth by bp ratio """
	sampled_bp = 0
	total_bp = 0
	count_reads = 0
	paths = mate1.split(',')
	if mate2: paths += mate2.split(',')
	divisor = 4 if file_type == 'fastq' else 2
	for path in paths:
		for line_number, line in enumerate(utility.iopen(path)):
			if (line_number-1) % divisor == 0:
				total_bp += len(line) - 1
				if count_reads < max_reads:
					sampled_bp += len(line) - 1
				count_reads += 1
	bp_ratio = float(total_bp)/sampled_bp
	for id in species_abundance:
		species_abundance[id]['cov'] *= bp_ratio
		species_abundance[id]['count'] = int(species_abundance[id]['count'] * bp_ratio)

def run_pipeline(args):
	
	""" Run entire pipeline """
	# read info files
	species_info = read_annotations(args)
	marker_info = read_marker_info(args)
		
	# align reads
	start = time()
	print("\nAligning reads to marker-genes database")
	args['log'].write("\nAligning reads to marker-genes database\n")
	map_reads_hsblast(args)
	print("  %s minutes" % round((time() - start)/60, 2))
	print("  %s Gb maximum memory" % utility.max_mem_usage())

	# find best hit for each read
	start = time()
	print("\nClassifying reads")
	args['log'].write("\nClassifying reads\n")
	best_hits = find_best_hits(args, marker_info)
	unique_alns = assign_unique(args, best_hits, species_info, marker_info)
	species_alns = assign_non_unique(args, best_hits, unique_alns, marker_info)
	print("  %s minutes" % round((time() - start)/60, 2))
	print("  %s Gb maximum memory" % utility.max_mem_usage())
	
	# estimate species abundance
	start = time()
	print("\nEstimating species abundance")
	args['log'].write("\nEstimating species abundance\n")
	total_gene_length = read_gene_lengths(args, species_info, marker_info)
	species_abundance = normalize_counts(species_alns, total_gene_length)
	print("  %s minutes" % round((time() - start)/60, 2) )
	print("  %s Gb maximum memory" % utility.max_mem_usage())
	
	# scale read depth
	if args['max_reads'] and args['scale']:
		start = time()
		print("\nScaling species read-depth to total metagenome size")
		args['log'].write("\nScaling species read-depth to total metagenome size\n")
		scale_depth(species_abundance, args['m1'], args['m2'], args['file_type'], args['max_reads'])
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())
	
	# write results
	write_abundance(args['outdir'], species_abundance, species_info)

	# clean up
	if args['remove_temp']:
		import shutil
		shutil.rmtree('%s/species/temp' % args['outdir'])
