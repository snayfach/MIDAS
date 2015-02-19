#!/usr/bin/python

# MicrobeCNV - estimation of gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3

__version__ = '0.0.1'

# TO DO
# import microbe_species module
# print cluster name, size, and relative abundance when read-mapping
# add option to 'pick-up' at different points in pipeline (use existing bam files)

# Libraries
# ---------
import sys
import os
import numpy as np
import argparse
import pysam
import gzip
import time
import subprocess
import operator
import Bio.SeqIO

# Functions
# ---------
def parse_arguments():
	""" Parse command line arguments """
	
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
	
	parser.add_argument('--version', action='version', version='MicrobeCNV %s' % __version__)
	parser.add_argument('-v', '--verbose', action='store_true', default=False)
	
	input = parser.add_argument_group('Input')
	input.add_argument('-1', type=str, dest='m1', help='FASTQ file containing 1st mate')
	input.add_argument('-2', type=str, dest='m2', help='FASTQ file containing 2nd mate')
	input.add_argument('-U', type=str, dest='r', help='FASTQ file containing unpaired reads')
	input.add_argument('-D', type=str, dest='db_dir', help='Directory of bt2 indexes for genome clusters')
	
	output = parser.add_argument_group('Output')
	output.add_argument('-o', type=str, dest='out', help='Base name for output files')

	pipe = parser.add_argument_group('Pipeline')
	pipe.add_argument('--all', action='store_true', dest='all',
		default=False, help='Run entire pipeline')
	pipe.add_argument('--profile', action='store_true', dest='profile',
		default=False, help='Fast estimation of genome-cluster abundance')
	pipe.add_argument('--align', action='store_true', dest='align',
		default=False, help='Align reads to genome-clusters')
	pipe.add_argument('--map', action='store_true', dest='map',
		default=False, help='Assign reads to mapping locations')
	pipe.add_argument('--cov', action='store_true', dest='cov',
		default=False, help='Compute coverage of pangenomes')
		
	aln = parser.add_argument_group('Alignment')
	aln.add_argument('--reads', type=int, dest='reads', help='Number of reads to use from sequence file (use all)')
	aln.add_argument('--abun', type=float, dest='abun', default=0.05,
			help='Abundance threshold for aligning to genome cluster (0.05)')
			
	map = parser.add_argument_group('Mapping')
	map.add_argument('--pid', type=float, dest='pid', default=90,
			help='Minimum percent identity between read and reference (90.0)')
	
	return vars(parser.parse_args())

def check_arguments(args):
	""" Check validity of command line arguments """
	
	# Pipeline options
	if not any([args['all'], args['profile'], args['align'], args['map'], args['cov']]):
		sys.exit('Specify pipeline option(s): --all, --profile, --align, --map, --cov')
	if args['all']:
		args['profile'] = True
		args['align'] = True
		args['map'] = True
		args['cov'] = True

	# Input options
	if (args['m1'] or args['m2']) and args['r']:
		sys.exit('Cannot use both -1/-2 and -U')
	if (args['m1'] and not args['m2']) or (args['m2'] and not args['m1']):
		sys.exit('Must specify both -1 and -2 for paired-end reads')
	if not (args['m1'] or args['r']):
		sys.exit('Specify reads using either -1 and -2 or -U')
	if args['m1'] and not os.path.isfile(args['m1']):
		sys.exit('Input file specified with -1 does not exist')
	if args['m2'] and not os.path.isfile(args['m2']):
		sys.exit('Input file specified with -2 does not exist')
	if args['r'] and not os.path.isfile(args['r']):
		sys.exit('Input file specified with -U does not exist')
	if args['db_dir'] and not os.path.isdir(args['db_dir']):
		sys.exit('Input directory specified with --db-dir does not exist')

	# Output options
	if not args['out']:
		sys.exit('Specify output basename with -o')

def print_copyright():
	# print out copyright information
	print ("-------------------------------------------------------------------------")
	print ("MicrobeCNV - estimation of gene-copy-number from shotgun sequence data")
	print ("version %s; github.com/snayfach/MicrobeCNV" % __version__)
	print ("Copyright (C) 2015 Stephen Nayfach")
	print ("Freely distributed under the GNU General Public License (GPLv3)")
	print ("-------------------------------------------------------------------------\n")

def parse_profile(inpath):
	""" Parse output from MicrobeSpecies """
	infile = open(inpath)
	next(infile)
	for line in infile:
		fields = [
			('cluster_id', str), ('mapped_reads', int), ('prop_mapped', float),
			('cell_count', float), ('prop_cells', float), ('avg_pid', float)]
		values = line.rstrip().split()
		yield dict( [ (f[0], f[1](v)) for f, v in zip(fields, values)] )

def select_genome_clusters(args):
	""" Select genome clusters to map to """
	cluster_to_abun = {}
	inpath = '%s.species' % args['out']
	if not os.path.isfile(inpath):
		sys.exit("Could not locate species profile: %s" % inpath)
	for rec in parse_profile(inpath):
		if rec['prop_cells'] >= args['abun']:
			cluster_to_abun[rec['cluster_id']] = rec['prop_cells']
	return cluster_to_abun

def align_reads(genome_clusters):
	""" Use Bowtie2 to map reads to all specified genome clusters """
	for cluster_id in genome_clusters:
		index_bn = '/'.join([args['db_dir'], cluster_id, cluster_id])
		# Build command
		command = 'bowtie2 --no-unal --very-sensitive -x %s ' % index_bn
		#   max reads to search
		if args['reads']: command += '-u %s ' % args['reads']
		#   input files
		if args['m1']: command += '-1 %s -2 %s ' % (args['m1'], args['m2'])
		else: command += '-U %(r)s ' % args['r']
		#   output
		command += '| samtools view -b - > %s' % '.'.join([args['out'], cluster_id, 'bam'])
		# Run command
		if args['verbose']: print("  running: %s") % command
		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = process.communicate()

def fetch_paired_reads(bam_path):
	""" Use pysam to yield paired end reads from bam file """
	pe_read = []
	aln_file = pysam.AlignmentFile(bam_path, "rb")
	for aln in aln_file.fetch(until_eof = True):
		if aln.mate_is_unmapped and aln.is_read1:
			yield [aln]
		elif aln.mate_is_unmapped and aln.is_read2:
			yield [aln]
		else:
			pe_read.append(aln)
			if len(pe_read) == 2:
				yield pe_read
				pe_read = []

def compute_aln_score(pe_read):
	""" Compute alignment score for paired-end read """
	if pe_read[0].mate_is_unmapped:
		score = pe_read[0].query_length - dict(pe_read[0].tags)['NM']
		return score
	else:
		score1 = pe_read[0].query_length - dict(pe_read[0].tags)['NM']
		score2 = pe_read[1].query_length - dict(pe_read[1].tags)['NM']
		return score1 + score2

def compute_perc_id(pe_read):
	""" Compute percent identity for paired-end read """
	if pe_read[0].mate_is_unmapped:
		length = pe_read[0].query_length
		edit = dict(pe_read[0].tags)['NM']
	else:
		length = pe_read[0].query_length + pe_read[1].query_length
		edit = dict(pe_read[0].tags)['NM'] + dict(pe_read[1].tags)['NM']
	return 100 * (length - edit)/float(length)

def find_best_hits(genome_clusters):
	""" Find top scoring alignment for each read """
	if args['verbose']: print("  finding best alignments across GCs...")
	best_hits = {}
	
	# map reads across genome clusters
	for cluster_id in genome_clusters:
		bam_path = '.'.join([args['out'], cluster_id, 'bam'])
		if not os.path.isfile(bam_path): # check that bam file exists
			sys.stderr.write("    bam file not found for genome-cluster %s. skipping\n" % cluster_id)
			continue
		if args['verbose']: print("    parsing: %s") % bam_path
		for pe_read in fetch_paired_reads(bam_path):
			# parse pe_read
			query = pe_read[0].query_name
			score = compute_aln_score(pe_read)
			pid = compute_perc_id(pe_read)
			if pid < args['pid']: # filter aln
				continue
			elif query not in best_hits: # store aln
				best_hits[query] = {'score':score, 'aln':{cluster_id:pe_read} }
			elif score > best_hits[query]['score']: # update aln
				best_hits[query] = {'score':score, 'aln':{cluster_id:pe_read} }
			elif score == best_hits[query]['score']: # append aln
				best_hits[query]['aln'][cluster_id] = pe_read
	return best_hits

def report_mapping_summary(best_hits):
	""" Summarize hits to genome-clusters """
	hit1, hit2, hit3 = 0, 0, 0
	for value in best_hits.values():
		if len(value['aln']) == 1: hit1 += 1
		elif len(value['aln']) == 2: hit2 += 1
		else: hit3 += 1
	print("  summary:")
	print("    %s reads assigned to any GC (%s)" % (hit1+hit2+hit3, round(float(hit1+hit2+hit3)/args['reads'], 2)) )
	print("    %s reads assigned to 1 GC (%s)" % (hit1, round(float(hit1)/args['reads'], 2)) )
	print("    %s reads assigned to 2 GCs (%s)" % (hit2, round(float(hit2)/args['reads'], 2)) )
	print("    %s reads assigned to 3 or more GCs (%s)" % (hit3, round(float(hit3)/args['reads'], 2)) )

	
def resolve_ties(best_hits, cluster_to_abun):
	""" Reassign reads that map equally well to >1 genome cluster """
	if args['verbose']: print("  reassigning reads mapped to >1 GC...")
	for query, rec in best_hits.items():
		if len(rec['aln']) == 1:
			best_hits[query] = rec['aln'].items()[0]
		if len(rec['aln']) > 1:
			target_gcs = rec['aln'].keys()
			abunds = [cluster_to_abun[gc] for gc in target_gcs]
			probs = [abund/sum(abunds) for abund in abunds]
			selected_gc = np.random.choice(target_gcs, 1, p=probs)[0]
			best_hits[query] = (selected_gc, rec['aln'][selected_gc])
	return best_hits

def write_best_hits(selected_clusters, best_hits, out):
	""" Write reassigned PE reads to disk """
	if args['verbose']: print("  writing mapped reads to disk...")
	# open filehandles
	aln_files = {}
	for cluster_id in selected_clusters:
		inpath = '.'.join([out, cluster_id, 'bam'])
		if not os.path.isfile(inpath): continue
		infile = pysam.AlignmentFile(inpath, 'rb')
		outpath = '.'.join([out, cluster_id, 'reassigned', 'bam'])
		aln_files[cluster_id] = pysam.AlignmentFile(outpath, 'wb', template=infile)
	# write reads to disk
	for cluster_id, pe_read in best_hits.values():
		for aln in pe_read:
			if aln is not None:
				aln_files[cluster_id].write(aln)

def compute_bed_cov(p_cov, p_bed):
    """ Compute coverage of features """
    # Read in coverage
    si_pos_to_cov = {}
    for line in open(p_cov):
        si, pos, cov = line.rstrip().split()
        si_pos_to_cov[si, int(pos)] = float(cov)
    # Compute coverage of features
    gene_to_cov = {}
    f_in = gzip.open(p_bed)
    for line in f_in:
        cov = 0
        si, start, stop, gene_oid = line.rstrip().split()
        for pos in range(int(start), int(stop)+1):
            cov += si_pos_to_cov[si, pos]
        gene_to_cov[gene_oid] = cov
    return gene_to_cov

# Main
# ------

args = parse_arguments()
check_arguments(args)

if args['verbose']: print_copyright()

if args['profile']:
	if args['verbose']: print("Estimating the abundance of genome-clusters")
	pass
cluster_to_abun = select_genome_clusters(args)
selected_clusters = cluster_to_abun.keys()

if args['align']:
	if args['verbose']: print("Aligning reads to reference genomes")
	align_reads(selected_clusters)

if args['map']:
	if args['verbose']: print("Mapping reads to most likely genomic positions")
	best_hits = find_best_hits(selected_clusters)
	if args['verbose']: report_mapping_summary(best_hits)
	best_hits = resolve_ties(best_hits, cluster_to_abun)
	write_best_hits(selected_clusters, best_hits, args['out'])

if args['cov']:
	if args['verbose']: print("Computing coverage of pangenomes")


#def aggregate_alignments(args, paths, alns):
#	""" Group all alignments to each genome cluster """
#	cluster_to_aln = dict([(x.rstrip(),[]) for x in open(paths['clusters']).readlines()])
#	if args['bam:
#		aln_file = pysam.AlignmentFile(args['bam if args['bam else '%s.bam' % args['out, "rb")
#		for aln in alns:
#			cluster_id, genome_id, gene_id, marker_id = aln_file.getrname(aln.reference_id).split('_')
#			cluster_to_aln[cluster_id].append(aln)
#	elif args['m8:
#		for aln in alns:
#			cluster_id, genome_id, gene_id, marker_id = aln['target'].split('_')
#			cluster_to_aln[cluster_id].append(aln)
#	return cluster_to_aln
#
#def compute_avg_mapq(args, alns):
#	""" Compute average map quality for list of alignments """
#	if len(alns) == 0:
#		return 'NA'
#	elif args['bam:
#		return np.mean([1 - dict(aln.tags)['NM']/float(aln.query_length) for aln in alns])
#	elif args['m8:
#		return 'NA'
#
#def compute_avg_pid(args, alns):
#	""" Compute average percent identity for list of alignments """
#	if len(alns) == 0:
#		return 'NA'
#	elif args['bam:
#		return np.mean([1 - dict(aln.tags)['NM']/float(aln.query_length) for aln in alns])
#	elif args['m8:
#		return np.mean([aln['pid']/float(100) for aln in alns])
#
#def alignment_summary(cluster_to_aln, args):
#	""" Write summary to outfile """
#	outfile = open(args['out+'.summary', 'w')
#	fields = ['cluster_id', 'mapped_reads', 'relabun', 'avg_pid', 'avg_mapq']
#	outfile.write('\t'.join(fields)+'\n')
#	total_mapped = sum([len(aln) for aln in cluster_to_aln.values()])
#	for cluster_id, alns in cluster_to_aln.items():
#		mapped_reads = len(alns)
#		relabun = len(alns)/float(total_mapped)
#		avg_pid = compute_avg_pid(args, alns)
#		avg_mapq = compute_avg_mapq(args, alns)
#		record = [str(x) for x in [cluster_id, mapped_reads, relabun, avg_pid, avg_mapq]]
#		outfile.write('\t'.join(record)+'\n')

