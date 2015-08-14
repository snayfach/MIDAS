#!/usr/bin/python

# Libraries
# ---------
import sys
import os
import gzip
import subprocess
import Bio.SeqIO
import random
import numpy as np
import time
from platform import system

# Functions
# ---------

def parse_relative_paths(args):
	""" Identify relative file and directory paths """
	paths = {}
	if system() not in ['Linux', 'Darwin']:
		sys.exit("Operating system '%s' not supported" % system())
	else:
		main_dir = os.path.dirname(os.path.abspath(__file__))
		paths['blastn'] = '/'.join([main_dir, 'bin', system(), 'blastn'])
		assert(os.path.isfile(paths['blastn']))
		paths['fa_to_fq'] = '/'.join([main_dir, 'bin', 'fa_to_fq.py'])
		assert(os.path.isfile(paths['blastn']))
		paths['cluster_ids'] = '/'.join([main_dir,'data','cluster_annotations.txt'])
		assert(os.path.isfile(paths['cluster_ids']))
		paths['gene_length'] = '/'.join([main_dir,'data','gene_length.txt'])
		assert(os.path.isfile(paths['gene_length']))
		paths['marker_cutoffs'] = '/'.join([main_dir,'data','pid_cutoffs.txt'])
		assert(os.path.isfile(paths['marker_cutoffs']))
		paths['db'] = '/'.join([main_dir,'data','phyeco.blastdb'])
	return paths

def quality_filter(seq, qual, args):
	""" Return true if read fails QC """
	l = len(seq)
	# check for Ns
	if sum([1 if b == 'N' else 0 for b in seq])/float(l) > args['max_n']:
		return True
	# check length
	if l < args['min_length']:
		return True
	# read passed QC
	return False

def map_reads_blast(args, paths):
	""" Use blastn to map reads in fasta file to marker database """
	command = 'python %s %s %s | ' % (paths['fa_to_fq'], args['inpath'], args['nreads']) # convert from fq to fq & pipe to blastn
	command += '%s ' % paths['blastn']
	command += '-query /dev/stdin -db %s ' % paths['db']
	command += '-out /dev/stdout -outfmt 6 -num_threads %s' % args['threads']
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = process.communicate()
	return stdout

def parse_blast(blastout):
	""" Yield formatted record from BLAST m8 file """
	formats = [str,str,float,int,float,float,float,float,float,float,float,float]
	fields = ['query','target','pid','aln','mis','gaps','qstart','qend','tstart','tend','evalue','score']
	for line in blastout.split('\n')[0:-1]:
		values = line.rstrip().split()
		yield dict([(field, format(value)) for field, format, value in zip(fields, formats, values)])

def find_best_hits(blastout, paths):
	""" Find top scoring alignment for each read """
	best_hits = {}
	marker_cutoffs = get_markers(paths)
	for aln in blastout:
		marker_id = aln['target'].split('_')[-1]
		if aln['pid'] < marker_cutoffs[marker_id]: # does not meet marker cutoff
			continue
		elif aln['query'] not in best_hits: # record aln
			best_hits[aln['query']] = [aln]
		elif best_hits[aln['query']][0]['score'] == aln['score']: # add aln
			best_hits[aln['query']] += [aln]
		elif best_hits[aln['query']][0]['score'] < aln['score']: # update aln
			best_hits[aln['query']] = [aln]
	return best_hits.values()

def assign_unique(args, paths, alns):
	""" Count the number of uniquely mapped reads to each genome cluster """
	unique_alns = dict([(x.rstrip().split()[0],[]) for x in open(paths['cluster_ids']).readlines()])
	unique = 0
	non_unique = 0
	for aln in alns:
		if len(aln) == 1:
			unique += 1
			cluster_id = aln[0]['target'].split('_')[0]
			unique_alns[cluster_id].append(aln[0])
		else:
			non_unique += 1
	if args['verbose']:
		print("\tuniquely mapped reads: %s" % unique)
		print("\tambiguously mapped reads: %s" % non_unique)
	return unique_alns

def assign_non_unique(args, paths, alns, unique_alns):
	""" Probabalistically assign ambiguously mapped reads """
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

def get_markers(paths):
	""" Read in optimal mapping parameters for marker genes """
	marker_cutoffs = {}
	infile = open(paths['marker_cutoffs'])
	for line in infile:
		marker_id, min_pid = line.rstrip().split()
		marker_cutoffs[marker_id] = float(min_pid)
	return marker_cutoffs

def impute_missing_args(args):
	""" Fill in defaults for missing arguments """
	if 'temp_dir' not in args:
		args['temp_dir'] = None
	if 'indir' not in args:
		args['indir'] = None
	if 'max_n' not in args:
		args['max_n'] = 1.0
	if 'min_quality' not in args:
		args['min_quality'] = 0
	if 'min_length' not in args:
		args['min_length'] = 0
	if 'verbose' not in args:
		args['verbose'] = False
	if 'threads' not in args:
		args['threads'] = 1
	if 'normalize' not in args:
		args['normalize'] = False
	if 'nreads' not in args:
		args['nreads'] = float('Inf')
	return args

def estimate_mix_props(alns, paths):
	""" Count the number of uniquely mapped reads to each genome cluster """
	unique_counts = dict([(x.rstrip().split()[0],0) for x in open(paths['cluster_ids']).readlines()])
	unique = 0
	non_unique = 0
	for aln in alns:
		if len(aln) == 1:
			unique += 1
			cluster_id = aln[0]['target'].split('_')[0]
			unique_counts[cluster_id] += 1
		else:
			non_unique += 1
	if args['verbose']:
		print("Uniquely mapped reads: %s" % unique)
		print("Ambiguously mapped reads: %s" % non_unique)
	return unique_counts

def read_gene_lengths(paths):
	""" Read in total gene length per cluster_id """
	total_gene_length = dict([(x.rstrip().split()[0],0) for x in open(paths['cluster_ids']).readlines()])
	for line in open(paths['gene_length']):
		cluster_id = line.split()[0].split('_')[0]
		gene_length = int(line.rstrip().split()[1])
		total_gene_length[cluster_id] += gene_length
	return total_gene_length

def normalize_counts(cluster_alns, total_gene_length):
	""" Normalize counts by gene length and sum contrain """
	# norm by gene length, compute rpkg, compute cov
	cluster_abundance = {}
	for cluster_id, alns in cluster_alns.items():
		cluster_abundance[cluster_id] = {}
		# compute coverage
		if len(alns) > 0:
			bp = sum([aln['aln'] for aln in alns])
			cov = float(bp)/total_gene_length[cluster_id]
		else:
			cov = 0.0
		# store results
		cluster_abundance[cluster_id] = {'cov':cov}
	# compute relative abundance
	total_cov = sum([_['cov'] for _ in cluster_abundance.values()])
	for cluster_id in cluster_abundance.keys():
		cov = cluster_abundance[cluster_id]['cov']
		cluster_abundance[cluster_id]['rel_abun'] = cov/total_cov if total_cov > 0 else 0
	return cluster_abundance

def write_abundance(outpath, cluster_abundance):
	""" Write cluster results to specified output file """
	outfile = open(outpath, 'w')
	fields = ['cluster_id', 'coverage', 'relative_abundance']
	outfile.write('\t'.join(fields)+'\n')
	for cluster_id, values in cluster_abundance.items():
		record = [cluster_id, values['cov'], values['rel_abun']]
		outfile.write('\t'.join([str(x) for x in record])+'\n')

def estimate_species_abundance(args):
	""" Run entire pipeline """

	# impute missing args & get relative file paths
	args = impute_missing_args(args)
	paths = parse_relative_paths(args)

	# align reads
	start = time.time()
	if args['verbose']: print("Aligning reads")
	blastout = parse_blast(map_reads_blast(args, paths))
	if args['verbose']: print("\t %ss" % round(time.time() - start))
	
	# find best hit for each read
	start = time.time()
	if args['verbose']: print("Classifying reads")
	best_hits = find_best_hits(blastout, paths)
	unique_alns = assign_unique(args, paths, best_hits)
	cluster_alns = assign_non_unique(args, paths, best_hits, unique_alns)
	if args['verbose']: print("\t %ss" % round(time.time() - start))
	
	# estimate genome cluster abundance
	start = time.time()
	if args['verbose']: print("Estimating cluster abundance")
	total_gene_length = read_gene_lengths(paths)
	cluster_abundance = normalize_counts(cluster_alns, total_gene_length)
	if args['verbose']: print("\t %ss" % round(time.time() - start))
	
	# return results
	return cluster_abundance



