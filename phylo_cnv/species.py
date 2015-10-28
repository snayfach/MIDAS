#!/usr/bin/python

# Libraries
# ---------
import sys, os, subprocess
import random
import numpy as np
from time import time
import platform
import tempfile
import misc

# Functions
# ---------

def add_tempfiles(args):
	""" Identify relative file and directory paths """
	args['tempfile'] = tempfile.mkstemp(dir=os.path.dirname(args['out']), suffix='.read_count')[1]
	args['blastout'] = tempfile.mkstemp(dir=os.path.dirname(args['out']), suffix='.m8')[1]


def map_reads_hsblast(args):
	""" Use hs-blastn to map reads in fasta file to marker database """
	# fasta to fastq
	command = 'python %s' % args['fa_to_fq']
	command += ' %s' % args['m1'] # fastq
	if args['m2']: command += ',%s' % args['m2'] # and mate if specified
	if args['reads']: command += ' %s' % args['reads'] # number of reads if specified
	command += ' 2> %s' % args['tempfile'] # tmpfile to store # of reads, bp sampled
	# hs-blastn
	command += ' | %s align' % args['hs-blastn']
	if args['speed'] == 'sensitive': command += ' -word_size 18' # decrease word size for more sensisitve search
	command += ' -query /dev/stdin -db %s/%s/hs-blast' % (args['db'], 'marker_genes') # specify db
	command += ' -outfmt 6 -num_threads %s' % args['threads'] # specify num threads
	command += ' -out %s' % args['blastout'] # output file
	command += ' -evalue 1e-3' # %id for reporting hits
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = process.communicate()

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
	for aln in parse_blast(args['blastout']):
		i += 1
		marker_id = aln['target'].split('_')[-1]
		if aln['pid'] < marker_cutoffs[marker_id]: # does not meet marker cutoff
			continue
		elif query_coverage(aln) < 0.75: # filter local alignments
			continue
		elif aln['query'] not in best_hits: # record aln
			best_hits[aln['query']] = [aln]
		elif best_hits[aln['query']][0]['score'] == aln['score']: # add aln
			best_hits[aln['query']] += [aln]
		elif best_hits[aln['query']][0]['score'] < aln['score']: # update aln
			best_hits[aln['query']] = [aln]
	if args['verbose']: print("  total alignments: %s" % i)
	return best_hits.values()

def assign_unique(args, alns):
	""" Count the number of uniquely mapped reads to each genome cluster """
	unique_alns = dict([(x.rstrip().split()[0],[]) for x in open(args['cluster_ids']).readlines()])
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
		print("  uniquely mapped reads: %s" % unique)
		print("  ambiguously mapped reads: %s" % non_unique)
	return unique_alns

def assign_non_unique(args, alns, unique_alns):
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

def get_markers(args):
	""" Read in optimal mapping parameters for marker genes """
	marker_cutoffs = {}
	infile = open(args['pid_cutoffs'])
	for line in infile:
		marker_id, min_pid = line.rstrip().split()
		marker_cutoffs[marker_id] = float(min_pid)
	return marker_cutoffs

def estimate_mix_props(alns, args):
	""" Count the number of uniquely mapped reads to each genome cluster """
	unique_counts = dict([(x.rstrip().split()[0],0) for x in open(args['cluster_ids']).readlines()])
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

def read_gene_lengths(args):
	""" Read in total gene length per cluster_id """
	total_gene_length = dict([(x.rstrip().split()[0],0) for x in open(args['cluster_ids']).readlines()])
	for line in open(args['gene_length']):
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

def estimate_abundance(args):
	
	""" Run entire pipeline """
	# impute missing args & get relative file paths
	add_tempfiles(args)
	misc.add_executables(args)
	misc.add_data_files(args)
	
	# align reads
	start = time()
	if args['verbose']: print("\nAligning reads to marker genes database")
	map_reads_hsblast(args)
	if args['verbose']:
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory") % misc.max_mem_usage()
	
	# find best hit for each read
	start = time()
	if args['verbose']: print("\nClassifying reads")
	best_hits = find_best_hits(args)
	unique_alns = assign_unique(args, best_hits)
	cluster_alns = assign_non_unique(args, best_hits, unique_alns)
	if args['verbose']:
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory") % misc.max_mem_usage()
	
	# estimate genome cluster abundance
	start = time()
	if args['verbose']: print("\nEstimating species abundance")
	total_gene_length = read_gene_lengths(args)
	cluster_abundance = normalize_counts(cluster_alns, total_gene_length)
	if args['verbose']:
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory") % misc.max_mem_usage()
	
	# convert to cellular relative abundances
	if args['norm']:
		start = time()
		if args['verbose']: print("\nConverting to cellular relative abundances")
		from microbe_census import microbe_census
		ags = microbe_census.run_pipeline({'seqfile':args['m1']})[0]
		reads, bp = [int(x) for x in open(args['tempfile']).read().rstrip().split()]
		genomes = bp/float(ags)
		for cluster_id in cluster_abundance:
			cov = cluster_abundance[cluster_id]['cov']
			cluster_abundance[cluster_id]['rel_abun'] = cov/genomes
		if args['verbose']:
			print("  average genome size: %s" % round(ags,2))
			print("  total bp sampled: %s" % bp)
			print("  total genome coverage: %s" % round(genomes,2))
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % misc.max_mem_usage()

	# write results
	write_abundance(args['out'], cluster_abundance)

	# clean up
	if not args['keep_temp']:
		os.remove(args['tempfile'])
		os.remove(args['blastout'])

def write_abundance(outpath, cluster_abundance):
	""" Write cluster results to specified output file """
	outfile = open(outpath, 'w')
	fields = ['cluster_id', 'coverage', 'relative_abundance']
	outfile.write('\t'.join(fields)+'\n')
	for cluster_id, values in cluster_abundance.items():
		record = [cluster_id, values['cov'], values['rel_abun']]
		outfile.write('\t'.join([str(x) for x in record])+'\n')

def read_abundance(inpath):
	""" Parse output from PhyloSpecies """
	if not os.path.isfile(inpath):
		sys.exit("Could not locate species profile: %s\nTry rerunning with --species_profile" % inpath)
	dict = {}
	fields = [('cluster_id', str), ('cov', float), ('rel_abun', float)]
	infile = open(inpath)
	next(infile)
	for line in infile:
		values = line.rstrip().split()
		dict[values[0]] = {}
		for field, value in zip(fields[1:], values[1:]):
			dict[values[0]][field[0]] = field[1](value)
	return dict

def select_genome_clusters(args):
	""" Select genome clusters to map to """
	import operator
	cluster_sets = {}
	# read in cluster abundance if necessary
	if any([args['gc_topn'], args['gc_cov'], args['gc_rbun']]):
		cluster_abundance = read_abundance(args['profile'])
		# user specifed a coverage threshold
		if args['gc_cov']:
			cluster_sets['gc_cov'] = set([])
			for cluster_id, values in cluster_abundance.items():
				if values['cov'] >= args['gc_cov']:
					cluster_sets['gc_cov'].add(cluster_id)
		# user specifed a relative-abundance threshold
		if args['gc_rbun']:
			cluster_sets['gc_rbun'] = set([])
			for cluster_id, values in cluster_abundance.items():
				if values['rel_abun'] >= args['gc_rbun']:
					cluster_sets['gc_rbun'].add(cluster_id)
		# user specifed topn genome-clusters
		if args['gc_topn']:
			cluster_sets['gc_topn'] = set([])
			cluster_abundance = [(i,d['rel_abun']) for i,d in cluster_abundance.items()]
			sorted_abundance = sorted(cluster_abundance, key=operator.itemgetter(1), reverse=True)
			for cluster_id, rel_abun in sorted_abundance[0:args['gc_topn']]:
				cluster_sets['gc_topn'].add(cluster_id)
	# user specified a list of one or more genome-clusters
	if args['gc_id']:
		cluster_sets['gc_rbun'] = set([])
		for cluster_id in args['gc_id']:
			cluster_sets['gc_rbun'].add(cluster_id)
	# intersect sets of genome-clusters
	my_clusters = list(set.intersection(*cluster_sets.values()))
	# check that specified genome-clusters are valid
	ref_db = os.path.join(args['db'], 'genome_clusters')
	for cluster_id in my_clusters:
		if cluster_id not in os.listdir(ref_db) :
			sys.exit("\nError: the specified genome_cluster '%s' was not found in the reference database:\n%s\n" % (cluster_id, ref_db))
	# remove bad cluster_ids
	for line in open(args['bad_gcs']):
		try: my_clusters.remove(line.rstrip())
		except: pass
	# check that at least one genome-cluster was selected
	if len(my_clusters) == 0:
		sys.exit("\nError: no genome-clusters sastisfied your selection criteria. \n")
	return my_clusters

