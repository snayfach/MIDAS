#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, subprocess
from time import time
from midas import utility

def read_annotations(args):
	annotations = {}
	inpath = '%s/species_info.txt' % args['db']
	if not os.path.isfile(inpath): sys.exit("File not found: %s" % inpath)
	for rec in utility.parse_file(inpath):
		annotations[rec['species_id']] = rec['species_name']
	return annotations

def map_reads_hsblast(args):
	""" Use hs-blastn to map reads in fasta file to marker database """
	# stream sequences
	command = 'python %s' % args['stream_seqs']
	command += ' -1 %s' % args['m1'] # fasta/fastq
	if args['m2']: command += ' -2 %s' % args['m2'] # mate
	if args['max_reads']: command += ' -n %s' % args['max_reads'] # number of reads
	if args['read_length']: command += ' -l %s' % args['read_length'] # read length
	command += ' 2> %s/species/temp/read_count.txt' % args['outdir'] # tmpfile to store # of reads, bp sampled
	# hs-blastn
	command += ' | %s align' % args['hs-blastn']
	command += ' -word_size %s' % args['word_size']
	command += ' -query /dev/stdin'
	command += ' -db %s/%s/%s' % (args['db'], 'marker_genes', args['db_type'])
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
	return float(aln['aln'])/int(aln['query'].split('_')[-1])

def find_best_hits(args):
	""" Find top scoring alignment for each read """
	best_hits = {}
	marker_cutoffs = get_markers(args)
	i = 0
	qcovs = []
	for aln in parse_blast('%s/species/temp/alignments.m8' % args['outdir']):
		i += 1
		marker_id = aln['target'].split('_')[-1]
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

def assign_unique(args, alns, species_ids):
	""" Count the number of uniquely mapped reads to each genome cluster """
	unique_alns = dict([(_,[]) for _ in species_ids])
	unique = 0
	non_unique = 0
	for aln in alns:
		if len(aln) == 1:
			unique += 1
			cluster_id = aln[0]['target'].split('_')[0]
			unique_alns[cluster_id].append(aln[0])
		else:
			non_unique += 1
	print("  uniquely mapped reads: %s" % unique)
	print("  ambiguously mapped reads: %s" % non_unique)
	return unique_alns

def assign_non_unique(args, alns, unique_alns):
	""" Probabalistically assign ambiguously mapped reads """
	import numpy as np
	import random
	total_alns = unique_alns.copy()
	for aln in alns:
		if len(aln) > 1:
			clusters = [x['target'].split('_')[0] for x in aln]
			counts = [len(unique_alns[x]) for x in clusters]
			if sum(counts) == 0:
				cluster_id = random.sample(clusters, 1)[0]
			else:
				probs = [float(count)/sum(counts) for count in counts]
				cluster_id = np.random.choice(clusters, 1, p=probs)[0]
			total_alns[cluster_id].append(aln[clusters.index(cluster_id)])
	return total_alns

def get_markers(args):
	""" Read in optimal mapping parameters for marker genes; override if user has provided cutoff """
	marker_cutoffs = {}
	inpath = '/'.join([args['db'], 'marker_genes/pid_cutoffs.txt'])
	if not os.path.isfile(inpath): sys.exit("File not found: %s" % inpath)
	for line in open(inpath):
		marker_id, min_pid = line.rstrip().split()
		if args['mapid']:
			marker_cutoffs[marker_id] = args['mapid']
		else:
			marker_cutoffs[marker_id] = float(min_pid)
	return marker_cutoffs

def read_gene_lengths(args, species_ids):
	""" Read in total gene length per cluster_id """
	total_gene_length = dict([(_,0) for _ in species_ids])
	inpath = '/'.join([args['db'], 'marker_genes/gene_length.txt'])
	if not os.path.isfile(inpath): sys.exit("File not found: %s" % inpath)
	for line in open(inpath):
		gene_id, gene_length, db_type = line.rstrip().split()
		cluster_id = gene_id.split('_')[0]
		if db_type == args['db_type']:
			total_gene_length[cluster_id] += int(gene_length)
	return total_gene_length

def normalize_counts(cluster_alns, total_gene_length):
	""" Normalize counts by gene length and sum contrain """
	# norm by gene length, compute rpkg, compute cov
	cluster_abundance = {}
	total_cov = 0.0
	for cluster_id, alns in cluster_alns.items():
		cluster_abundance[cluster_id] = {}
		# compute coverage
		if len(alns) > 0:
			bp = sum([aln['aln'] for aln in alns])
			cov = float(bp)/total_gene_length[cluster_id]
		else:
			cov = 0.0
		# store results
		cluster_abundance[cluster_id] = {'cov':cov, 'count':len(alns)}
		total_cov += cov
	# compute relative abundance
	total_cov = sum([_['cov'] for _ in cluster_abundance.values()])
	for cluster_id in cluster_abundance.keys():
		cov = cluster_abundance[cluster_id]['cov']
		cluster_abundance[cluster_id]['rel_abun'] = cov/total_cov if total_cov > 0 else 0
	print("  total marker-gene coverage: %s" % round(total_cov, 3))
	return cluster_abundance

def estimate_abundance(args):
	
	""" Run entire pipeline """
	# impute missing args & get relative file paths
	species_info = read_annotations(args)
	
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
	best_hits = find_best_hits(args)
	unique_alns = assign_unique(args, best_hits, species_info)
	cluster_alns = assign_non_unique(args, best_hits, unique_alns)
	print("  %s minutes" % round((time() - start)/60, 2))
	print("  %s Gb maximum memory" % utility.max_mem_usage())
	
	# estimate genome cluster abundance
	start = time()
	print("\nEstimating species abundance")
	args['log'].write("\nEstimating species abundance\n")
	total_gene_length = read_gene_lengths(args, species_info)
	species_abundance = normalize_counts(cluster_alns, total_gene_length)
	print("  %s minutes" % round((time() - start)/60, 2) )
	print("  %s Gb maximum memory" % utility.max_mem_usage())
	
	# write results
	write_abundance(args['outdir'], species_abundance, species_info)

	# clean up
	if args['remove_temp']:
		import shutil
		shutil.rmtree('%s/species/temp' % args['outdir'])

def write_abundance(outdir, cluster_abundance, annotations):
	""" Write cluster results to specified output file """
	outpath = '%s/species/species_profile.txt' % outdir
	outfile = open(outpath, 'w')
	fields = ['species_id', 'species_name', 'count_reads', 'coverage', 'relative_abundance']
	outfile.write('\t'.join(fields)+'\n')
	for cluster_id, values in cluster_abundance.items():
		record = [cluster_id, annotations[cluster_id], values['count'], values['cov'], values['rel_abun']]
		outfile.write('\t'.join([str(x) for x in record])+'\n')

def read_abundance(inpath):
	""" Parse species abundance file """
	if not os.path.isfile(inpath):
		sys.exit("\nCould not locate species profile: %s\nTry rerunning with run_midas.py species" % inpath)
	abun = {}
	for rec in utility.parse_file(inpath):
		# format record
		if 'cluster_id' in rec: rec['species_id'] = rec['cluster_id']
		if 'count_reads' in rec: rec['count_reads'] = int(rec['count_reads'])
		if 'coverage' in rec: rec['coverage'] = float(rec['coverage'])
		if 'relative_abundance' in rec: rec['relative_abundance'] = float(rec['relative_abundance'])
		abun[rec['species_id']] = rec
	return abun

def select_genome_clusters(args):
	""" Select genome clusters to map to """
	import operator
	cluster_sets = {}
	# read in cluster abundance if necessary
	if any([args['gc_topn'], args['gc_cov']]):
		cluster_abundance = read_abundance('%s/species/species_profile.txt' % args['outdir'])
		# user specifed a coverage threshold
		if args['gc_cov']:
			cluster_sets['gc_cov'] = set([])
			for cluster_id, values in cluster_abundance.items():
				if values['coverage'] >= args['gc_cov']:
					cluster_sets['gc_cov'].add(cluster_id)
		# user specifed topn genome-clusters
		if args['gc_topn']:
			cluster_sets['gc_topn'] = set([])
			cluster_abundance = [(i,d['relative_abundance']) for i,d in cluster_abundance.items()]
			sorted_abundance = sorted(cluster_abundance, key=operator.itemgetter(1), reverse=True)
			for cluster_id, rel_abun in sorted_abundance[0:args['gc_topn']]:
				cluster_sets['gc_topn'].add(cluster_id)
	# user specified a list of one or more genome-clusters
	if args['gc_id']:
		cluster_sets['gc_id'] = set([])
		for cluster_id in args['gc_id']:
			cluster_sets['gc_id'].add(cluster_id)
	# intersect sets of genome-clusters
	my_clusters = list(set.intersection(*cluster_sets.values()))
	# check that specified genome-clusters are valid
	ref_db = os.path.join(args['db'], 'genome_clusters')
	for cluster_id in my_clusters:
		if cluster_id not in os.listdir(ref_db) :
			sys.exit("\nError: the specified species_id '%s' was not found in the reference database:\n%s\n" % (cluster_id, ref_db))
	# remove bad cluster_ids
	inpath = '/'.join([args['db'], 'exclude.txt'])
	if not os.path.isfile(inpath): sys.exit("File not found: %s" % inpath)
	for line in open(inpath):
		try: my_clusters.remove(line.rstrip())
		except: pass
	# check that at least one genome-cluster was selected
	if len(my_clusters) == 0:
		sys.exit("\nError: no species sastisfied your selection criteria. \n")
	return my_clusters

