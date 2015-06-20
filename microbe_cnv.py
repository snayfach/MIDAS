#!/usr/bin/python

# MicrobeCNV - estimation of gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3

__version__ = '0.0.1'

# TO DO
# integrate microbe_species
# compute %id of pangenes

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
from math import ceil

# Functions
# ---------
def parse_arguments():
	""" Parse command line arguments """
	
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
	
	parser.add_argument('--version', action='version', version='MicrobeCNV %s' % __version__)
	parser.add_argument('-v', '--verbose', action='store_true', default=False)
	
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-1', type=str, dest='m1', help='FASTQ file containing 1st mate if paired or unpaired reads')
	io.add_argument('-2', type=str, dest='m2', help='FASTQ file containing 2nd mate if paired')
	io.add_argument('-D', type=str, dest='db_dir', help='Directory of bt2 indexes for genome clusters')
	io.add_argument('-o', type=str, dest='out', help='Directory for output files')

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
	pipe.add_argument('--extract', action='store_true', dest='extract',
		default=False, help='Extract mapped reads from bam file & write to FASTQ')
	pipe.add_argument('--remap', action='store_true', dest='remap',
		default=False, help='Re-map reads to representative genomes')
	pipe.add_argument('--snps', action='store_true', dest='snps',
		default=False, help='Run samtools mpileup & estimate SNP frequencies')
				
	gc = parser.add_argument_group('Genome-cluster inclusion (choose one)')
	gc.add_argument('--gc_topn', type=int, dest='gc_topn', default=5, help='Top N most abundant (5)')
	gc.add_argument('--gc_cov', type=float, dest='gc_cov', help='Coverage threshold (None)')
	gc.add_argument('--gc_rbun', type=float, dest='gc_rbun', help='Relative abundance threshold (None)')
	gc.add_argument('--gc_id', type=str, dest='gc_id', help='Identifier of specific genome cluster (None)')
	gc.add_argument('--gc_list', type=str, dest='gc_list', help='Comma-separated list of genome cluster ids (None)')
	
	reads = parser.add_argument_group('Read selection')
	reads.add_argument('--rd_ms', type=int, dest='reads_ms',
		default=5000000, help='# reads to use for estimating genome-cluster abundance (5,000,000)')
	reads.add_argument('--rd_align', type=int, dest='reads_align',
		help='# reads to use for pangenome alignment (All)')
	reads.add_argument('--rd_batch', type=int, dest='rd_batch', default=5000000,
		help='Batch size in # reads. Smaller batch sizes requires less memory, but can take longer to run (5,000,000)')

	map = parser.add_argument_group('Mapping')
	map.add_argument('--map_pid', type=float, dest='pid',
		default=90, help='Minimum percent identity between read and reference (90.0)')
		
	snps = parser.add_argument_group('SNP detection')
	snps.add_argument('--snps_mapq', type=str, dest='snps_mapq',
		default='0', help='Minimum map quality (0)')
	snps.add_argument('--snps_baseq', type=str, dest='snps_baseq',
		default='0', help='Minimum base quality (0)')
				
	mask = parser.add_argument_group('Leave-One-Out (for simulated data only)')
	mask.add_argument('--tax_mask', action='store_true', dest='tax_mask',
		default=False, help='Discard alignments for reads and ref seqs from the same genome')
	mask.add_argument('--tax_map', type=str, dest='tax_map',
		help='File mapping read ids to genome ids')
	
	return vars(parser.parse_args())

def check_arguments(args):
	""" Check validity of command line arguments """
	
	# Pipeline options
	if not any([args['all'], args['profile'], args['align'], args['map'],
	            args['cov'], args['extract'], args['remap'], args['snps']]):
		sys.exit('Specify one or more pipeline option(s): --all, --profile, --align, --map, --cov, --extract, --remap, --snps')
	if args['all']:
		args['profile'] = True
		args['align'] = True
		args['map'] = True
		args['cov'] = True
		args['extract'] = True
		args['remap'] = True
		args['snps'] = True
	if args['tax_mask'] and not args['tax_map']:
		sys.exit('Specify file mapping read ids in FASTQ file to genome ids in reference database')

	# Input options
	if not args['m1']:
		sys.exit('Specify input FASTQ file(s) with -1 -2 or -U')
	if args['m1'] and not os.path.isfile(args['m1']):
		sys.exit('Input file specified with -1 does not exist')
	if args['m2'] and not os.path.isfile(args['m2']):
		sys.exit('Input file specified with -2 does not exist')
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
	print ("-------------------------------------------------------------------------")

def read_microbe_species(inpath):
	""" Parse output from MicrobeSpecies """
	if not os.path.isfile(inpath):
		sys.exit("Could not locate species profile: %s\nTry rerunning with --profile" % inpath)
	dict = {}
	fields = [
		('cluster_id', str), ('reads', float), ('bp', float), ('rpkg', float),
		('cov', float), ('prop_cov', float), ('rel_abun', float)]
	infile = open(inpath)
	next(infile)
	for line in infile:
		values = line.rstrip().split()
		dict[values[0]] = {}
		for field, value in zip(fields[1:], values[1:]):
			dict[values[0]][field[0]] = field[1](value)
	return dict

def select_genome_clusters(cluster_abundance, args):
	""" Select genome clusters to map to """
	my_clusters = {}
	# user specified a single genome-cluster
	if args['gc_id']:
		cluster_id = args['gc_id']
		if cluster_id not in cluster_abundance:
			sys.exit("Error: specified genome-cluster id %s not found" % cluster_id)
		else:
			abundance = cluster_abundance[args['gc_id']]['rel_abun']
			my_clusters[args['gc_id']] = abundance
	# user specified a list of genome-clusters
	elif args['gc_list']:
		for cluster_id in args['gc_list'].split(','):
			if cluster_id not in cluster_abundance:
				sys.exit("Error: specified genome-cluster id %s not found" % cluster_id)
			else:
				abundance = cluster_abundance[cluster_id]['rel_abun']
				my_clusters[cluster_id] = coverage
	# user specifed a coverage threshold
	elif args['gc_cov']:
		for cluster_id, values in cluster_abundance.items():
			if values['cell_count'] >= args['gc_cov']:
				my_clusters[cluster_id] = values['cov']
	# user specifed a relative-abundance threshold
	elif args['gc_rbun']:
		for cluster_id, values in cluster_abundance.items():
			if values['prop_mapped'] >= args['gc_rbun']:
				my_clusters[cluster_id] = values['rel_abun']
	# user specifed a relative-abundance threshold
	elif args['gc_topn']:
		cluster_abundance = [(i,d['rel_abun']) for i,d in cluster_abundance.items()]
		sorted_abundance = sorted(cluster_abundance, key=operator.itemgetter(1), reverse=True)
		for cluster_id, coverage in sorted_abundance[0:args['gc_topn']]:
			my_clusters[cluster_id] = coverage
	return my_clusters

def align_reads(genome_clusters, batch_index, reads_start, batch_size, tax_mask):
	""" Use Bowtie2 to map reads to all specified genome clusters """
	# Create output directory
	outdir = os.path.join(args['out'], 'bam')
	if not os.path.isdir(outdir): os.mkdir(outdir)
	for cluster_id in genome_clusters:
		# Build command
		command = '%s --no-unal --very-sensitive ' % args['bowtie2']
		#   index
		command += '-x %s ' % '/'.join([args['db_dir'], cluster_id, 'btdb', cluster_id])
		#   specify reads
		command += '-s %s -u %s ' % (reads_start, batch_size)
		#	report up to 20 hits/read if masking hits
		if tax_mask: command += '-k 20 '
		#   input
		if (args['m1'] and args['m2']): command += '-1 %s -2 %s ' % (args['m1'], args['m2'])
		else: command += '-U %s' % args['m1']
		#   output
		bampath = '/'.join([args['out'], 'bam', '%s.%s.bam' % (cluster_id, batch_index)])
		command += '| %s view -b - > %s' % (args['samtools'], bampath)
		# Run command
		if args['verbose']: print("    running: %s") % command
		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = process.communicate()

def align_to_rep(genome_clusters):
	""" Use Bowtie2 to map reads to representative genomes from each genome cluster
	"""
	# Create output directory
	outdir = os.path.join(args['out'], 'bam_rep')
	if not os.path.isdir(outdir): os.mkdir(outdir)
	for cluster_id in genome_clusters:
		# Build command
		#	bowtie2
		command = '%s --no-unal --very-sensitive ' % args['bowtie2']
		#   bt2 index
		command += '-x %s ' % '/'.join([args['db_dir'], cluster_id, 'btdb_rep', cluster_id])
		#   input fastq
		command += '-U %s ' % '/'.join([args['out'], 'fastq', '%s.fastq.gz' % cluster_id])
		#   convert to bam
		command += '| %s view -b - ' % args['samtools']
		#   sort bam
		command += '| %s sort -f - %s ' % (args['samtools'], '/'.join([outdir, '%s.bam' % cluster_id]))
		# Run command
		if args['verbose']: print("    running: %s") % command
		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = process.communicate()

def pileup_on_rep(genome_clusters):
	""" Use Samtools to create pileup, filter low quality bases, and write results to VCF file """
	# Create output directory
	outdir = os.path.join(args['out'], 'vcf')
	if not os.path.isdir(outdir): os.mkdir(outdir)
	for cluster_id in genome_clusters:
		# Build command
		#   mpileup
		command = '%s mpileup -uv -A -d 10000 --skip-indels -B ' % args['samtools']
		#   quality filtering
		command += '-q %s -Q %s ' % (args['snps_mapq'], args['snps_baseq'])
		#   reference fna file
		command += '-f %s ' % '/'.join([args['db_dir'], cluster_id, 'representative.fna'])
		#   input bam file
		command += '%s ' % '/'.join([args['out'], 'bam_rep', '%s.bam' % cluster_id])
		#   output vcf file
		command += '> %s ' % '/'.join([outdir, '%s.vcf' % cluster_id])
		# Run command
		if args['verbose']: print("    running: %s") % command
		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = process.communicate()

def fetch_reads(aln_file):
	""" Use pysam to yield paired end reads from bam file """
	pe_read = []
	for aln in aln_file.fetch(until_eof = True):
		if not aln.is_paired:
			yield [aln]
		elif aln.mate_is_unmapped and aln.is_read1:
			yield [aln]
		elif aln.mate_is_unmapped and aln.is_read2:
			yield [aln]
		else:
			pe_read.append(aln)
			if len(pe_read) == 2:
				yield pe_read
				pe_read = []

def compute_aln_score(pe_read):
	""" Compute alignment score for single or paired-end read """
	if not pe_read[0].is_paired:
		score = pe_read[0].query_length - dict(pe_read[0].tags)['NM']
		return score
	elif pe_read[0].mate_is_unmapped:
		score = pe_read[0].query_length - dict(pe_read[0].tags)['NM']
		return score
	else:
		score1 = pe_read[0].query_length - dict(pe_read[0].tags)['NM']
		score2 = pe_read[1].query_length - dict(pe_read[1].tags)['NM']
		return score1 + score2

def compute_perc_id(pe_read):
	""" Compute percent identity for paired-end read """
	if not pe_read[0].is_paired:
		length = pe_read[0].query_length
		edit = dict(pe_read[0].tags)['NM']
	elif pe_read[0].mate_is_unmapped:
		length = pe_read[0].query_length
		edit = dict(pe_read[0].tags)['NM']
	else:
		length = pe_read[0].query_length + pe_read[1].query_length
		edit = dict(pe_read[0].tags)['NM'] + dict(pe_read[1].tags)['NM']
	return 100 * (length - edit)/float(length)

def find_best_hits(genome_clusters, batch_index, tax_mask, tax_map):
	""" Find top scoring alignment(s) for each read """
	if args['verbose']: print("    finding best alignments across GCs")
	best_hits = {}
	reference_map = {} # (cluster_id, ref_index) = ref_id (ref_id == scaffold id)
	
	# map reads across genome clusters
	for cluster_id in genome_clusters:
	
		# if masking alignments, read in:
		if tax_mask:
			scaffold_to_genome = {} # 1) map of scaffold to genome id
			inpath = '/'.join([args['db_dir'], cluster_id, 'genome_to_scaffold.gz'])
			infile = gzip.open(inpath)
			for line in infile:
				genome_id, scaffold_id = line.rstrip().split('\t')
				scaffold_to_genome[scaffold_id] = genome_id
			run_to_genome = {} # 2) run_accession to genome id
			for line in open(tax_map):
				run_accession, genome_id = line.rstrip().split('\t')
				run_to_genome[run_accession] = genome_id
		
		# get path to bam file
		bam_path = '/'.join([args['out'], 'bam', '%s.%s.bam' % (cluster_id, batch_index)])
		if not os.path.isfile(bam_path):
			sys.stderr.write("      warning: bam file not found for %s.%s" % (cluster_id, batch_index))
			continue
			
		# loop over PE reads
		aln_file = pysam.AlignmentFile(bam_path, "rb")
		for pe_read in fetch_reads(aln_file):
		
			# map reference ids
			for aln in pe_read:
				ref_index = aln.reference_id
				ref_id = aln_file.getrname(ref_index).split('|')[1] # reformat ref id
				reference_map[(cluster_id, ref_index)] = ref_id
				
			# mask alignment
			if tax_mask:
				ref_index = pe_read[0].reference_id
				ref_id = aln_file.getrname(ref_index).split('|')[1]
				run_accession = pe_read[0].query_name.split('.')[0]
				if run_to_genome[run_accession] == scaffold_to_genome[ref_id]:
					continue
					
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
				
	# resolve ties
	best_hits = resolve_ties(best_hits, genome_clusters)
	return best_hits, reference_map

def report_mapping_summary(best_hits):
	""" Summarize hits to genome-clusters """
	hit1, hit2, hit3 = 0, 0, 0
	for value in best_hits.values():
		if len(value['aln']) == 1: hit1 += 1
		elif len(value['aln']) == 2: hit2 += 1
		else: hit3 += 1
	if args['reads_align']:
		print("  summary:")
		print("    %s reads assigned to any GC (%s)" % (hit1+hit2+hit3, round(float(hit1+hit2+hit3)/args['reads_align'], 2)) )
		print("    %s reads assigned to 1 GC (%s)" % (hit1, round(float(hit1)/args['reads_align'], 2)) )
		print("    %s reads assigned to 2 GCs (%s)" % (hit2, round(float(hit2)/args['reads_align'], 2)) )
		print("    %s reads assigned to 3 or more GCs (%s)" % (hit3, round(float(hit3)/args['reads_align'], 2)) )
	else:
		print("  summary:")
		print("    %s reads assigned to any GC" % (hit1+hit2+hit3))
		print("    %s reads assigned to 1 GC" % (hit1))
		print("    %s reads assigned to 2 GCs" % (hit2))
		print("    %s reads assigned to 3 or more GCs" % (hit3))

def resolve_ties(best_hits, cluster_to_abun):
	""" Reassign reads that map equally well to >1 genome cluster """
	if args['verbose']: print("    reassigning reads mapped to >1 GC")
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

def write_best_hits(genome_clusters, best_hits, reference_map, batch_index):
	""" Write reassigned PE reads to disk """
	if args['verbose']: print("    writing mapped reads to disk")
	try: os.makedirs('/'.join([args['out'], 'reassigned']))
	except: pass
	# open filehandles
	aln_files = {}
	scaffold_to_genome = {}
	# loop over genome clusters
	for cluster_id in genome_clusters:
		# get template bam file
		bam_path = '/'.join([args['out'], 'bam', '%s.%s.bam' % (cluster_id, batch_index)])
		if not os.path.isfile(bam_path):
			sys.stderr.write("    bam file not found for %s.%s Skipping\n" % (cluster_id, batch_index))
			continue
		template = pysam.AlignmentFile(bam_path, 'rb')
		# store filehandle
		outpath = '/'.join([args['out'], 'reassigned', '%s.%s.bam' % (cluster_id, batch_index)])
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

def compute_pangenome_coverage(cluster_id, batch_index, read_length):
	""" Use bedtools to compute coverage of pangenome """
	global pangene_to_cov
	bedcov_out = run_bed_coverage(cluster_id, batch_index) # run bedtools
	for r in parse_bed_cov(bedcov_out): # aggregate coverage by pangene_id
		pangene_id = r['pangene_id']
		coverage = r['reads'] * read_length / r['gene_length']
		pangene_to_cov[pangene_id] += coverage

def run_bed_coverage(cluster_id, batch_index):
	""" Run bedCoverage for cluster_id """
	bampath = '/'.join([args['out'], 'reassigned', '%s.%s.bam' % (cluster_id, batch_index)])
	bedpath = '/'.join([args['db_dir'], cluster_id, '%s.bed' % cluster_id])
	cmdargs = {'bedcov':args['bedcov'], 'bam':bampath, 'bed':bedpath}
	command = '%(bedcov)s -abam %(bam)s -b %(bed)s' % cmdargs
	process = subprocess.Popen(command % cmdargs, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	return out

def compute_phyeco_cov(pangene_to_cov, cluster_id):
	""" Compute coverage of phyeco markers for genome cluster """
	phyeco_covs = []
	inpath = '/'.join([args['db_dir'], cluster_id, 'pangene_to_phyeco.gz'])
	infile = gzip.open(inpath)
	next(infile)
	for line in infile:
		pangene, phyeco_id = line.rstrip().split()
		phyeco_covs.append(pangene_to_cov[pangene])
	return np.median(phyeco_covs)

def max_mem_usage():
	""" Return max mem usage (Gb) of self and child processes """
	max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
	return round((max_mem_self + max_mem_child)/float(1e6), 2)

def get_read_length(args):
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
	return np.mean(read_lengths)

def get_read_count(inpath):
	""" Count the number of reads in fastq file """
	line_count = 0
	infile = gzip.open(inpath)
	for line_count, line in enumerate(infile):
		pass
	return (line_count+1)/4

def batch_reads(args, batch_size):
	""" Define batches of reads (batch_index, reads_start, reads_end)"""
	batches = []
	total_reads = args['reads_align'] if args['reads_align'] else get_read_count(args['m1'])
	nbatches = int(ceil(total_reads/float(batch_size)))
	for batch_index in range(1, nbatches+1):
		reads_start = (batch_index * batch_size) - batch_size + 1
		reads_end = reads_start + batch_size - 1
		if reads_end > total_reads: # adjust size for final chunk
			batch_size = batch_size - (reads_end - total_reads)
		batches.append([batch_index, reads_start, batch_size])
	return batches

def convert_to_ascii_quality(scores):
	""" Convert quality scores to Sanger encoded (Phred+33) ascii values """
	ascii = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQR"""
	score_to_ascii = dict((x,y) for x,y in zip(range(0,50),list(ascii)))
	return ''.join([score_to_ascii[x] for x in scores])

def convert_from_ascii_quality(asciis):
	""" Convert quality scores to Sanger encoded (Phred+33) ascii values """
	ascii = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQR"""
	ascii_to_score = dict((y,x) for x,y in zip(range(0,50),list(ascii)))
	return [ascii_to_score[x] for x in asciis]

def write_fastq_record(aln, index, outfile):
	""" Write pysam alignment record to outfile in FASTQ format """
	outfile.write('@%s.%s length=%s\n' % (aln.query_name,str(index),str(aln.query_length)))
	outfile.write('%s\n' % (aln.query_sequence))
	outfile.write('+%s.%s length=%s\n' % (aln.query_name,str(index),str(aln.query_length)))
	outfile.write('%s\n' % convert_to_ascii_quality(aln.query_qualities))

def pileup_to_bases(pileup, ref_base):
	""" Takes pileup string from mpileup output. Returns the aligned bases """
	bases = []
	index = 0
	while True:
		# parse pileup
		pileup_char = pileup[index]
		if pileup_char == '^': # skip '^' (start of a read segment) and ASCII character following '^'
			index += 2
		elif pileup_char == '$': # skip '$' (end of a read segment)
			index += 1
		elif pileup_char in ['+', '-']: # skip insertions and deletions (max 999 bp)
			for ndigits in [3,2,1]:
				indel_length = pileup[index+1:index+1+ndigits]
				if indel_length.isdigit(): break
			index += 1 + ndigits + int(indel_length)
		else: # keep positions for aligned bases
			if pileup_char in ['.',',']:
				bases.append(ref_base.upper())
			else:
				bases.append(pileup_char.upper())
			index += 1
		# stop when end of pileup reached
		if index == len(pileup):
			break
	return bases

def parse_pileup(inpath):
	""" Yields formatted records from mpileup output file """
	infile = gzip.open(inpath)
	for line in infile:
		r = line.rstrip().split()
		ref_id = r[0]
		ref_pos = int(r[1])
		ref_base = r[2].upper()
		depth = int(r[3])
		if depth == 0:
			bases = []
			phred_qualities = []
			map_qualities = []
		else:
			bases = pileup_to_bases(r[4], ref_base)
			phred_qualities = convert_from_ascii_quality(r[5])
			map_qualities = convert_from_ascii_quality(r[6])
		yield {'ref_id':ref_id, 'ref_pos':ref_pos,
		       'ref_base':ref_base, 'depth':depth,
			   'bases':bases, 'phred_qualities':phred_qualities,
			   'map_qualities':map_qualities}

def bam_to_fastq(genome_clusters, args):
	""" Converts bam to fastq for reads assigned to each genome cluster """
	bam_dir = '/'.join([args['out'], 'reassigned'])
	batch_indexes = set([_.split('.')[1] for _ in os.listdir(bam_dir)])
	fastq_dir = '/'.join([args['out'], 'fastq'])
	try: os.mkdir(fastq_dir)
	except: pass
	for genome_cluster in genome_clusters:
		outfile = gzip.open(os.path.join(fastq_dir, genome_cluster+'.fastq.gz'), 'w')
		for batch_index in batch_indexes:
			bam_name = '.'.join([genome_cluster, batch_index, 'bam'])
			bam_path = '/'.join([bam_dir, bam_name])
			aln_file = pysam.AlignmentFile(bam_path, "rb")
			for index, aln in enumerate(aln_file.fetch(until_eof = True)):
				write_fastq_record(aln, index, outfile)

def pileup_to_snps(genome_clusters, args):
	""" Parse pileups in order to call consensus alleles and reference allele frequencies """
	# create outdir
	outdir = '/'.join([args['out'], 'snps'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	for cluster_id in genome_clusters:
		# open outfile
		outfile = gzip.open('/'.join([args['out'], 'snps', '%s.snps.gz' % cluster_id]), 'w')
		header = ['ref_id', 'ref_pos', 'ref_base', 'coverage', 'counts',
				  'cons_base', 'cons_freq', 'mean_mapq', 'mean_baseq']
		outfile.write('\t'.join(header)+'\n')
		# parse pileup
		inpath = '/'.join([args['out'], 'pileup', '%s.pileup.gz' % cluster_id])
		for r in parse_pileup(inpath):
			if r['depth'] == 0:
				counts = [0, 0, 0, 0]
				cons_freq, cons_base, mean_mapq, mean_baseq = 'NA', 'NA', 'NA', 'NA'
			else:
				counts = [r['bases'].count(x) for x in ['A','T','C','G']]
				cons_freq = max(counts)/float(r['depth'])
				cons_base = ['A','T','C','G'][counts.index(max(counts))]
				mean_mapq = np.mean(r['map_qualities'])
				mean_baseq = np.mean(r['phred_qualities'])
			rec = [r['ref_id'], r['ref_pos'], r['ref_base'], r['depth'],
				   ','.join([str(x) for x in counts]),
				   cons_base, cons_freq, mean_mapq, mean_baseq]
			outfile.write('\t'.join([str(x) for x in rec])+'\n')


def vcf_to_snps(genome_clusters):
	""" Parse vcf file in order to call consensus alleles and reference allele frequencies """
	# create outdir
	outdir = '/'.join([args['out'], 'snps'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	for cluster_id in genome_clusters:
		# open outfile
		outfile = gzip.open('/'.join([args['out'], 'snps', '%s.snps.gz' % cluster_id]), 'w')
		header = ['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
		          'count_alleles', 'count_ref', 'count_alt', 'depth', 'ref_freq']
		outfile.write('\t'.join(header)+'\n')
		# parse vcf
		inpath = '/'.join([args['out'], 'vcf', '%s.vcf' % cluster_id])
		for r in parse_vcf(inpath):
			rec = [r['ref_id'], r['ref_pos'], r['ref_allele'], r['alt_allele'], r['cons_allele'],
				   r['count_alleles'], r['count_ref'], r['count_alt'], r['depth'], r['ref_freq']]
			outfile.write('\t'.join([str(x) for x in rec])+'\n')

def parse_vcf(inpath):
	""" Yields formatted records from VCF output file """
	infile = open(inpath)
	for line in infile:
		# skip header and split line
		if line[0] == '#': continue
		r = line.rstrip().split()
		# get alt alleles
		alt_alleles = r[4].split(',')
		if '<X>' in alt_alleles: alt_alleles.remove('<X>')
		count_alleles = 1 + len(alt_alleles)
		# get allele counts
		info = dict([(_.split('=')) for _ in r[7].split(';')])
		counts = [int(_) for _ in info['I16'].split(',')[0:4]]
		# get consensus allele
		# *note: occassionally there are counts for alternate alleles, but no listed alternate alleles
		if sum(counts) == 0:
			cons_allele = 'NA'
		elif sum(counts[0:2]) >= sum(counts[2:4]):
			cons_allele = r[3]
		elif len(alt_alleles) == 0:
			cons_allele = 'NA'
		else:
			cons_allele = alt_alleles[0]
		# yield formatted record
		yield {'ref_id':r[0],
			   'ref_pos':r[1],
			   'ref_allele':r[3],
			   'count_alleles':count_alleles,
			   'alt_allele':alt_alleles[0] if count_alleles > 1 else 'NA',
			   'depth':sum(counts),
			   'count_ref':sum(counts[0:2]),
			   'count_alt':sum(counts[2:4]),
			   'cons_allele':cons_allele,
			   'ref_freq':'NA'if sum(counts) == 0 else sum(counts[0:2])/float(sum(counts))
			   }





# Main
# ------

if __name__ == "__main__":

	args = parse_arguments()
	check_arguments(args) # need to check gc args
	
	src_dir = os.path.dirname(os.path.abspath(__file__))
	args['bowtie2'] = '/'.join([src_dir, 'lib', 'bowtie2-2.2.4', 'bowtie2'])
	args['samtools'] = '/'.join([src_dir, 'lib', 'samtools-1.2', 'samtools'])
	args['bcftools'] = '/'.join([src_dir, 'lib', 'bcftools-1.2', 'bcftools'])
	args['bedcov'] = '/'.join([src_dir, 'lib', 'bedtools2', 'bin', 'coverageBed'])
	
	if args['verbose']: print_copyright()

	if args['profile']:
		start = time.time()
		if args['verbose']: print("\nEstimating the abundance of genome-clusters")
		cluster_abundance = microbe_species.estimate_species_abundance(
			{'inpaths':[args['m1']], 'nreads':args['reads_ms'], 'outpath':'/'.join([args['out'], 'genome_clusters']),
			 'min_quality': 30, 'min_length': 50, 'max_n':0.05})
		microbe_species.write_results(
			'/'.join([args['out'], 'genome_clusters.abundance']),
			cluster_abundance
			)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	if args['verbose']: print("\nSelecting genome-clusters for pangenome alignment")
	cluster_abundance = read_microbe_species('/'.join([args['out'], 'genome_clusters.abundance']))
	genome_clusters = select_genome_clusters(cluster_abundance, args)
	if len(genome_clusters) == 0:
		sys.exit("No genome-clusters were detected that exceeded the minimum abundance threshold of %s" % args['abun'])
	elif args['verbose']:
		for cluster, abundance in sorted(genome_clusters.items(), key=operator.itemgetter(1), reverse=True):
			print("  cluster_id: %s abundance: %s" % (cluster, round(abundance,2)))

	if args['align']:
		start = time.time()
		if args['verbose']: print("\nAligning reads to reference genomes")
		for batch_index, reads_start, batch_size in batch_reads(args, args['rd_batch']):
			if args['verbose']: print("  batch %s:" % batch_index)
			align_reads(genome_clusters, batch_index, reads_start, batch_size, args['tax_mask'])
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	if args['map']:
		start = time.time()
		if args['verbose']: print("\nMapping reads to genome clusters")
		# get batch indexes from bam directory
		bam_dir = '/'.join([args['out'], 'bam'])
		batch_indexes = set([_.split('.')[1] for _ in os.listdir(bam_dir)])
		# loop over batch indexes
		for batch_index in batch_indexes:
			if args['verbose']: print("  batch %s:" % batch_index)
			best_hits, reference_map = find_best_hits(genome_clusters, batch_index, args['tax_mask'], args['tax_map'])
			write_best_hits(genome_clusters, best_hits, reference_map, batch_index)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	if args['cov']:
		start = time.time()
		if args['verbose']: print("\nEstimating coverage of pangenomes")
		# estimate average read length
		read_length = get_read_length(args)
		# get batch indexes from reassigned directory
		mapped_dir = '/'.join([args['out'], 'reassigned'])
		batch_indexes = set([_.split('.')[1] for _ in os.listdir(mapped_dir)])
		# loop over genome-clusters
		for cluster_id in genome_clusters:
			pangene_to_cov = defaultdict(float)
			if args['verbose']: print("  genome-cluster %s" % cluster_id)
			for batch_index in batch_indexes:
				compute_pangenome_coverage(cluster_id, batch_index, read_length)
			phyeco_cov = compute_phyeco_cov(pangene_to_cov, cluster_id)
			write_pangene_coverage(pangene_to_cov, phyeco_cov, cluster_id)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	if args['extract']:
		start = time.time()
		if args['verbose']: print("\nExtacting/writing mapped reads to FASTQ")
		bam_to_fastq(genome_clusters, args)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	if args['remap']:
		if args['verbose']: print("\nMapping reads to rep-genomes")
		start = time.time()
		align_to_rep(genome_clusters)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	if args['snps']:
		start = time.time()
		if args['verbose']: print("\nEstimating allele frequencies from mpileup")
		pileup_on_rep(genome_clusters)
		vcf_to_snps(genome_clusters)
		if args['verbose']:
			print("  %s minutes" % round((time.time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

#	if args['remap']:
#		start = time.time()
#		if args['verbose']: print("\nExtacting/writing mapped reads to FASTQ")
#		bam_to_fastq(genome_clusters, args)
#		if args['verbose']:
#			print("  %s minutes" % round((time.time() - start)/60, 2) )
#			print("  %s Gb maximum memory") % max_mem_usage()
#
#		if args['verbose']: print("\nMapping reads to rep-genomes & running mpileup")
#		start = time.time()
#		pileup_on_rep_genomes(genome_clusters)
#		if args['verbose']:
#			print("  %s minutes" % round((time.time() - start)/60, 2) )
#			print("  %s Gb maximum memory") % max_mem_usage()
#
#	if args['snps']:
#		start = time.time()
#		if args['verbose']: print("\nEstimating allele frequencies from mpileup")
#		pileup_to_snps(genome_clusters, args)
#		if args['verbose']:
#			print("  %s minutes" % round((time.time() - start)/60, 2) )
#			print("  %s Gb maximum memory") % max_mem_usage()
