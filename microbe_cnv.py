#!/usr/bin/python

# MicrobeCNV - estimation of gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3

__version__ = '0.0.1'

# TO DO
# time each section of code
# compare counting speed to samtools
# import microbe_species module
# print cluster name, size, and relative abundance when read-mapping
# compute coverage, %id of pangenes

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
import microbe_species
import resource
from collections import defaultdict

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
	output.add_argument('-o', type=str, dest='out', help='Directory for output files')

	pipe = parser.add_argument_group('Pipeline')
	pipe.add_argument('--all', action='store_true', dest='all',
		default=False, help='Run entire pipeline')
	pipe.add_argument('--profile', action='store_true', dest='profile',
		default=False, help='Estimate genome-cluster abundance using MicrobeSpecies')
	pipe.add_argument('--align', action='store_true', dest='align',
		default=False, help='Align reads to genome-clusters')
	pipe.add_argument('--map', action='store_true', dest='map',
		default=False, help='Assign reads to mapping locations')
	pipe.add_argument('--cov', action='store_true', dest='cov',
		default=False, help='Compute coverage of pangenomes')
		
	ms = parser.add_argument_group('Profiling Genome-Clusters')
	ms.add_argument('--reads_ms', type=int, dest='reads_ms',
		default=5000000, help='Number of reads to use for MicrobeSpecies (5,000,000)')
	
	aln = parser.add_argument_group('Alignment')
	aln.add_argument('--reads_bt', type=int, dest='reads_bt',
		help='Number of reads to use for genome alignment (use all)')
	aln.add_argument('--abun', type=float, dest='abun',
		default=1, help='Coverage threshold for aligning to genome cluster (1)')
			
	map = parser.add_argument_group('Mapping')
	map.add_argument('--pid', type=float, dest='pid',
		default=90, help='Minimum percent identity between read and reference (90.0)')
	
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
		sys.exit('Specify output directory with -o')

def print_copyright():
	# print out copyright information
	print ("-------------------------------------------------------------------------")
	print ("MicrobeCNV - estimation of gene-copy-number from shotgun sequence data")
	print ("version %s; github.com/snayfach/MicrobeCNV" % __version__)
	print ("Copyright (C) 2015 Stephen Nayfach")
	print ("Freely distributed under the GNU General Public License (GPLv3)")
	print ("-------------------------------------------------------------------------\n")

def read_microbe_species(inpath):
	""" Parse output from MicrobeSpecies """
	if not os.path.isfile(inpath):
		sys.exit("Could not locate species profile: %s" % inpath)
	dict = {}
	fields = [
		('cluster_id', str), ('mapped_reads', int), ('prop_mapped', float),
		('cell_count', float), ('prop_cells', float), ('avg_pid', float)]
	infile = open(inpath)
	next(infile)
	for line in infile:
		values = line.rstrip().split()
		dict[values[0]] = {}
		for field, value in zip(fields[1:], values[1:]):
			dict[values[0]][field[0]] = field[1](value)
	return dict

def select_genome_clusters(cluster_abundance):
	""" Select genome clusters to map to """
	my_clusters = {}
	for cluster_id, values in cluster_abundance.items():
		if values['cell_count'] >= args['abun']:
			my_clusters[cluster_id] = values['cell_count']
	return my_clusters

def align_reads(genome_clusters):
	""" Use Bowtie2 to map reads to all specified genome clusters """
	# Create output directory
	try: os.mkdir(os.path.join(args['out'], 'bam'))
	except: pass
	for cluster_id in genome_clusters:
		# Build command
		index_bn = '/'.join([args['db_dir'], cluster_id, cluster_id])
		command = '%s --no-unal --very-sensitive -x %s ' % (args['bowtie2'], index_bn)
		#   max reads to search
		if args['reads_bt']: command += '-u %s ' % args['reads_bt']
		#   input files
		if args['m1']: command += '-1 %s -2 %s ' % (args['m1'], args['m2'])
		else: command += '-U %(r)s ' % args['r']
		#   output
		command += '| %s view -b - > %s' % (args['samtools'], '/'.join([args['out'], 'bam', '%s.bam' % cluster_id]))
		# Run command
		if args['verbose']: print("  running: %s") % command
		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = process.communicate()

def fetch_paired_reads(aln_file):
	""" Use pysam to yield paired end reads from bam file """
	pe_read = []
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
	if args['verbose']: print("  finding best alignments across GCs:")
	best_hits = {}
	reference_map = {} # (cluster_id, ref_index) = ref_id (i think ref_id is the scaffold id)
	# map reads across genome clusters
	for cluster_id in genome_clusters:
		bam_path = '/'.join([args['out'], 'bam', '%s.bam' % cluster_id])
		if not os.path.isfile(bam_path): # check that bam file exists
			sys.stderr.write("    bam file not found for genome-cluster %s. skipping\n" % cluster_id)
			continue
		if args['verbose']: print("     %s") % os.path.basename(bam_path)
		aln_file = pysam.AlignmentFile(bam_path, "rb")
		for pe_read in fetch_paired_reads(aln_file):
			# map reference ids
			for aln in pe_read:
				ref_index = aln.reference_id
				ref_id = aln_file.getrname(ref_index).split('|')[1]
				reference_map[(cluster_id, ref_index)] = ref_id
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
	return best_hits, reference_map

def report_mapping_summary(best_hits):
	""" Summarize hits to genome-clusters """
	hit1, hit2, hit3 = 0, 0, 0
	for value in best_hits.values():
		if len(value['aln']) == 1: hit1 += 1
		elif len(value['aln']) == 2: hit2 += 1
		else: hit3 += 1
	if args['reads_bt']:
		print("  summary:")
		print("    %s reads assigned to any GC (%s)" % (hit1+hit2+hit3, round(float(hit1+hit2+hit3)/args['reads_bt'], 2)) )
		print("    %s reads assigned to 1 GC (%s)" % (hit1, round(float(hit1)/args['reads_bt'], 2)) )
		print("    %s reads assigned to 2 GCs (%s)" % (hit2, round(float(hit2)/args['reads_bt'], 2)) )
		print("    %s reads assigned to 3 or more GCs (%s)" % (hit3, round(float(hit3)/args['reads_bt'], 2)) )
	else:
		print("  summary:")
		print("    %s reads assigned to any GC" % (hit1+hit2+hit3))
		print("    %s reads assigned to 1 GC" % (hit1))
		print("    %s reads assigned to 2 GCs" % (hit2))
		print("    %s reads assigned to 3 or more GCs" % (hit3))

def resolve_ties(best_hits, cluster_to_abun):
	""" Reassign reads that map equally well to >1 genome cluster """
	if args['verbose']: print("  reassigning reads mapped to >1 GC")
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

def write_best_hits(selected_clusters, best_hits, reference_map):
	""" Write reassigned PE reads to disk """
	if args['verbose']: print("  writing mapped reads to disk")
	try: os.makedirs('/'.join([args['out'], 'reassigned']))
	except: pass
	# open filehandles
	aln_files = {}
	scaffold_to_genome = {}
	for bam_file in os.listdir('/'.join([args['out'], 'bam'])):
		# template bamfile
		cluster_id = bam_file.split('.')[0]
		inpath = '/'.join([args['out'], 'bam', bam_file])
		template = pysam.AlignmentFile(inpath, 'rb')
		# store filehandle
		outpath = '/'.join([args['out'], 'reassigned', '%s.bam' % cluster_id])
		aln_files[cluster_id] = pysam.AlignmentFile(outpath, 'wb', template=template)
	# write reads to disk
	for cluster_id, pe_read in best_hits.values():
		for aln in pe_read:
			aln_files[cluster_id].write(aln)

def write_pangene_coverage(pangene_to_cov, phyeco_cov, cluster_id):
	""" Write coverage of pangenes for genome cluster to disk """
	outdir = '/'.join([args['out'], 'coverage'])
	try: os.mkdir(outdir)
	except: pass
	outfile = gzip.open('/'.join([outdir, '%s.cov.gz' % cluster_id]), 'w')
	for pangene in sorted(pangene_to_cov.keys()):
		cov = pangene_to_cov[pangene]
		cn = cov/phyeco_cov if phyeco_cov > 0 else 0
		outfile.write('\t'.join([pangene, str(cov), str(cn)])+'\n')

def parse_bed_cov(bedcov_out):
	""" Yield dictionary of formatted values from bed coverage output """
	fields =  ['sid', 'start', 'end', 'gene_id', 'pangene_id', 'reads', 'pos_cov', 'gene_length', 'fract_cov']
	formats = [str, int, int, str, str, int, int, int, float]
	for line in bedcov_out.rstrip().split('\n'):
		rec = line.split()
		yield dict([(fields[i],formats[i](j)) for i,j in enumerate(rec)])

def compute_pangenome_coverage(cluster_id, read_length):
	""" Use bedtools to compute coverage of pangenome """
	bedcov_out = run_bed_coverage(cluster_id) # run bedtools
	pangene_to_cov = defaultdict(float) # init dict
	for r in parse_bed_cov(bedcov_out): # aggregate coverage by pangene_id
		pangene_id = r['pangene_id']
		coverage = r['reads'] * read_length / r['gene_length']
		pangene_to_cov[pangene_id] += coverage
	return pangene_to_cov

def run_bed_coverage(cluster_id):
	""" Run bedCoverage for cluster_id """
	bampath = '/'.join([args['out'], 'bam', '%s.bam' % cluster_id])
	bedpath = '/'.join([args['db_dir'], cluster_id, '%s.bed' % cluster_id])
	cmdargs = {'bedcov':args['bedcov'], 'bam':bampath, 'bed':bedpath}
	command = '%(bedcov)s -abam %(bam)s -b %(bed)s' % cmdargs
	process = subprocess.Popen(command % cmdargs, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	return out

def compute_phyeco_cov(pangene_to_cov, cluster_id):
	""" Compute coverage of phyeco markers for genome cluster """
	phyeco_covs = []
	inpath = '/'.join([args['db_dir'], cluster_id, '%s.phyeco.gz' % cluster_id])
	infile = gzip.open(inpath)
	next(infile)
	for line in infile:
		pangene_id, type, phyeco_id = line.rstrip().split()
		pangene = '_'.join([pangene_id, type])
		phyeco_covs.append(pangene_to_cov[pangene])
	return np.median(phyeco_covs)

def max_mem_usage():
	""" Return max mem usage (Gb) of self and child processes """
	max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
	return round((max_mem_self + max_mem_child)/float(1e6), 2)

def get_read_length():
	""" Estimate the average read length of fastq file from bam file """
	max_reads = 50000
	read_lengths = []
	bam_dir = '/'.join([args['out'], 'bam'])
	for file in os.listdir(bam_dir):
		bam_path = '/'.join([bam_dir, file])
		aln_file = pysam.AlignmentFile(bam_path, "rb")
		for index, aln in enumerate(aln_file.fetch(until_eof = True)):
			if index == max_reads: break
			else: read_lengths.append(aln.query_length)
	return np.median(read_lengths)


# Main
# ------

if __name__ == "__main__":

	args = parse_arguments()
	check_arguments(args)
	
	src_dir = os.path.dirname(os.path.abspath(__file__))
	args['bowtie2'] = '/'.join([src_dir, 'lib', 'bowtie2-2.2.4', 'bowtie2'])
	args['samtools'] = '/'.join([src_dir, 'lib', 'samtools-1.1', 'samtools'])
	args['bedcov'] = '/'.join([src_dir, 'lib', 'bedtools2', 'bin', 'coverageBed'])
	
	if args['verbose']: print_copyright()

	if args['profile']:
		start = time.time()
		if args['verbose']: print("Estimating the abundance of genome-clusters")
		cluster_abundance = microbe_species.estimate_species_abundance(
			{'inpaths': [args['m1']], 'nreads': args['reads_ms']})
		microbe_species.write_results(os.path.join(args['out'], 'cluster_abundance.txt'), cluster_abundance)
		selected_clusters = select_genome_clusters(cluster_abundance)
		if args['verbose']:
			for cluster, abundance in sorted(selected_clusters.items(), key=operator.itemgetter(1), reverse=True):
				print("  cluster_id: %s abundance: %s" % (cluster, round(abundance,2)))
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory\n") % max_mem_usage()
	else:
		cluster_abundance = read_microbe_species(os.path.join(args['out'], 'cluster_abundance.txt'))
		selected_clusters = select_genome_clusters(cluster_abundance)
	if len(selected_clusters) == 0:
		sys.exit("No genome-clusters were detected that exceeded the minimum abundance threshold of %s" % args['abun'])

	if args['align']:
		start = time.time()
		if args['verbose']: print("Aligning reads to reference genomes")
		align_reads(selected_clusters.keys())
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory\n") % max_mem_usage()

	if args['map']:
		start = time.time()
		if args['verbose']: print("Mapping reads to genome clusters")
		best_hits, reference_map = find_best_hits(selected_clusters.keys())
		if args['verbose']: report_mapping_summary(best_hits)
		best_hits = resolve_ties(best_hits, selected_clusters)
		write_best_hits(selected_clusters.keys(), best_hits, reference_map)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory\n") % max_mem_usage()

	if args['cov']:
		start = time.time()
		if args['verbose']: print("Computing coverage of pangenomes")
		read_length = get_read_length()
		for file in os.listdir('/'.join([args['out'], 'reassigned'])):
			cluster_id = file.split('.')[0]
			if args['verbose']: print("  %s") % cluster_id
			pangene_to_cov = compute_pangenome_coverage(cluster_id, read_length)
			phyeco_cov = compute_phyeco_cov(pangene_to_cov, cluster_id)
			write_pangene_coverage(pangene_to_cov, phyeco_cov, cluster_id)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory\n") % max_mem_usage()


