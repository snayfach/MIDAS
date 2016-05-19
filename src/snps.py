#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

# Libraries
# ---------
import sys
import os
import subprocess
import gzip
from time import time
import platform
import utility

# Functions
# ---------

def genome_align(args):
	""" Use Bowtie2 to map reads to representative genomes from each genome cluster
	"""
	# Build command
	#	bowtie2
	command = '%s --no-unal ' % args['bowtie2']
	#   index
	command += '-x %s ' % '/'.join([args['outdir'], 'snps/temp/genomes'])
	#   specify reads
	if args['max_reads']: command += '-u %s ' % args['max_reads']
	#   trim reads
	if args['trim']: command += '--trim3 %s ' % args['trim']
	#   speed/sensitivity
	command += '--%s ' % args['speed']
	#   threads
	command += '--threads %s ' % args['threads']
	#   file type
	if args['file_type'] == 'fasta': command += '-f '
	else: command += '-q '
	#   input file
	if (args['m1'] and args['m2']): command += '-1 %s -2 %s ' % (args['m1'], args['m2'])
	else: command += '-U %s' % args['m1']
	#   convert to bam
	command += '| %s view -b - ' % args['samtools']
	#   sort bam
	bam_path = os.path.join(args['outdir'], 'snps/temp/genomes.bam')
	command += '| %s sort -f - %s ' % (args['samtools'], bam_path)
	# Run command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# Check for errors
	utility.check_exit_code(process, command)
	utility.check_bamfile(args, bam_path)

def pileup(args):
	""" Filter alignments by % id, use samtools to create pileup, filter low quality bases, and write results to VCF file """
	# Build command
	#   percent id filtering
	bam_path = os.path.join(args['outdir'], 'snps/temp/genomes.bam')
	command  = 'python %s %s %s %s %s | ' % (args['stream_bam'], bam_path, '/dev/stdout', args['mapid'], args['readq'])
	#   mpileup
	command += '%s mpileup -uv -A -d 10000 --skip-indels ' % args['samtools']
	#   quality filtering
	if not args['baq']: command += '-B '
	#   quality filtering
	if args['redo_baq']: command += '-E '
	#   adjust MQ
	if args['adjust_mq']: command += '-C 50 '
	#   quality filtering
	command += '-q %s -Q %s ' % (args['mapq'], args['baseq'])
	#   reference fna file
	command += '-f %s ' % ('%s/snps/temp/genomes.fa' % args['outdir'])
	#   input bam file
	command += '- '
	#   output vcf file
	command += '> %s ' % ('%s/snps/temp/genomes.vcf' % args['outdir'])
	# Run command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)

def split_vcf(args):
	""" Format vcf output for easy parsing """
	inpath = os.path.join(args['outdir'], 'snps/temp/genomes.map')
	ref_to_cluster = utility.read_ref_to_cluster(inpath)
	# open outfiles for each cluster_id
	outdir = '/'.join([args['outdir'], 'snps/temp/vcf'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	outfiles = {}
	for cluster_id in set(ref_to_cluster.values()):
		outfiles[cluster_id] = open('/'.join([outdir, '%s.vcf' % cluster_id]), 'w')
	# parse vcf into temorary vcf files for each cluster_id
	for line in open('/'.join([args['outdir'], 'snps/temp/genomes.vcf'])):
		if line[0] == '#': continue
		cluster_id = ref_to_cluster[line.split()[0]]
		outfiles[cluster_id].write(line)
	# close outfiles
	for file in outfiles.values():
		file.close()

def read_ref_bases(args, cluster_id):
	""" Read in reference genome by position """
	import Bio.SeqIO
	ref = []
	centroid_path = '/'.join([args['db'], 'genome_clusters', cluster_id, 'representative.fna.gz'])
	infile = gzip.open(centroid_path)
	for rec in Bio.SeqIO.parse(infile, 'fasta'):
		for pos in range(1, len(rec.seq)+1):
			ref.append([rec.id, pos, rec.seq[pos-1].upper()])
	return sorted(ref)

def write_snp_header(outfile):
	""" Write header for formatted SNP file """
	fields = ['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
			  'count_alleles', 'count_ref', 'count_alt', 'depth', 'ref_freq']
	outfile.write('\t'.join(fields)+'\n')

def write_snp_record(outfile, snp, ref):
	""" Write record for formatted SNP file """
	if ref:
		snp = {'ref_id': ref[0], 'ref_pos': str(ref[1]), 'ref_allele': ref[2],
		       'alt_allele': 'NA', 'cons_allele': 'NA', 'count_alleles': 1,
		       'depth': 0, 'count_ref': 0, 'count_alt': 0, 'ref_freq': 'NA'}
	fields = ['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
			  'count_alleles', 'count_ref', 'count_alt', 'depth', 'ref_freq']
	record = [str(snp[field]) for field in fields]
	outfile.write('\t'.join(record)+'\n')

def format_vcf(args):
	""" Format vcf files to snp files and fill in missing positions """
	inpath = os.path.join(args['outdir'], 'snps/temp/genomes.map')
	ref_to_cluster = utility.read_ref_to_cluster(inpath)
	for cluster_id in set(ref_to_cluster.values()):
		# open outfile
		outfile = gzip.open('/'.join([args['outdir'], 'snps/output/%s.snps.gz' % cluster_id]), 'w')
		write_snp_header(outfile)
		# read sorted reference
		ref = read_ref_bases(args, cluster_id)
		ref_index = 0
		ref_length = len(ref)
		# write formatted records
		vcf_path = '/'.join([args['outdir'], 'snps/temp/vcf/%s.vcf' % cluster_id])
		for snp in parse_vcf(vcf_path): # loop over formatted records from vcf
			snp_pos = [snp['ref_id'], int(snp['ref_pos'])]
			while snp_pos != ref[ref_index][0:2]: # fill in missing snp positions
				write_snp_record(outfile, None, ref[ref_index]) # write missing record
				ref_index += 1
			write_snp_record(outfile, snp, None) # write present record
			ref_index += 1
		while ref_index < ref_length: # fill in trailing snps
			write_snp_record(outfile, None, ref[ref_index]) # write trailing record
			ref_index += 1

def parse_vcf(inpath):
	""" Yields formatted records from VCF output """
	infile = open(inpath)
	for line in infile:
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
			   'ref_freq':'NA' if sum(counts) == 0 else sum(counts[0:2])/float(sum(counts))
			   }

def snps_summary(args):
	""" Get summary of mapping statistics """
	# store stats
	stats = {}
	inpath = os.path.join(args['outdir'], 'snps/temp/genomes.map')
	ref_to_cluster = utility.read_ref_to_cluster(inpath)
	for cluster_id in set(ref_to_cluster.values()):
		genome_length, covered_bases, total_depth, identity, maf = [0,0,0,0,0]
		for r in utility.parse_file('/'.join([args['outdir'], 'snps/output/%s.snps.gz' % cluster_id])):
			genome_length += 1
			depth = int(r['depth'])
			if depth > 0:
				covered_bases += 1
				total_depth += depth
				if r['ref_allele'] == r['cons_allele']:
					identity += 1
				ref_freq = float(r['ref_freq'])
				maf += ref_freq if ref_freq <= 0.5 else 1 - ref_freq
		stats[cluster_id] = {'genome_length':genome_length,
							 'fraction_covered':covered_bases/float(genome_length),
							 'average_depth':total_depth/float(covered_bases) if covered_bases > 0 else 0,
							 'average_identity':identity/float(covered_bases) if covered_bases > 0 else 0,
							 'average_maf':maf/float(covered_bases) if covered_bases > 0 else 0
							 }
	# write stats
	fields = ['genome_length', 'fraction_covered', 'average_depth', 'average_identity', 'average_maf']
	outfile = open('/'.join([args['outdir'], 'snps/summary.txt']), 'w')
	outfile.write('\t'.join(['cluster_id'] + fields)+'\n')
	for cluster_id in stats:
		record = [cluster_id] + [str(stats[cluster_id][field]) for field in fields]
		outfile.write('\t'.join(record)+'\n')

def fetch_centroid(args, cluster_id):
	""" Get the genome_id corresponding to cluster centroid """
	inpath = '/'.join([args['db'], 'genome_clusters', cluster_id, 'genomes.txt.gz'])
	infile = gzip.open(inpath)
	for line in infile:
		if line.split()[2] == 'Y':
			return line.split()[1]

def build_genome_db(args, genome_clusters):
	""" Build FASTA and BT2 database from genome cluster centroids """
	# fasta database
	genomes_fasta = open('/'.join([args['outdir'], 'snps/temp/genomes.fa']), 'w')
	genomes_map = open('/'.join([args['outdir'], 'snps/temp/genomes.map']), 'w')
	db_stats = {'total_length':0, 'total_seqs':0, 'genome_clusters':0}
	for cluster_id in genome_clusters:
		if args['tax_mask'] and fetch_centroid(args, cluster_id) in args['tax_mask']:
			continue
		db_stats['genome_clusters'] += 1
		inpath = '/'.join([args['db'], 'genome_clusters', cluster_id, 'representative.fna.gz'])
		infile = gzip.open(inpath)
		for line in infile:
			genomes_fasta.write(line)
			db_stats['total_length'] += len(line.rstrip())
			if line[0] == '>':
				sid = line.rstrip().lstrip('>').split()[0]
				genomes_map.write(sid+'\t'+cluster_id+'\n')
				db_stats['total_seqs'] += 1
	# print out database stats
	print("  total genomes: %s" % db_stats['genome_clusters'])
	print("  total contigs: %s" % db_stats['total_seqs'])
	print("  total base-pairs: %s" % db_stats['total_length'])
	# bowtie2 database
	inpath = '/'.join([args['outdir'], 'snps/temp/genomes.fa'])
	outpath = '/'.join([args['outdir'], 'snps/temp/genomes'])
	command = ' '.join([args['bowtie2-build'], inpath, outpath])
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)

def remove_tmp(args):
	""" Remove specified temporary files """
	import shutil
	shutil.rmtree('/'.join([args['outdir'], 'snps/temp']))

def run_pipeline(args):
	""" Run entire pipeline """
	
	utility.add_executables(args) # Add paths to external files and binaries
	utility.add_data_files(args)
	
	# Build genome database for selected GCs
	if args['build_db']:
		import species
		print("\nBuilding database of representative genomes")
		start = time()
		genome_clusters = species.select_genome_clusters(args)
		build_genome_db(args, genome_clusters)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory") % utility.max_mem_usage()

	# Use bowtie2 to map reads to a representative genome for each genome-cluster
	if args['align']:
		args['file_type'] = utility.auto_detect_file_type(args['m1'])
		print("\nMapping reads to representative genomes")
		start = time()
		genome_align(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory") % utility.max_mem_usage()

	# Use mpileup to identify SNPs
	if args['call']:
		start = time()
		print("\nRunning mpileup")
		pileup(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory") % utility.max_mem_usage()

	# Split vcf into files for each GC, format, and report summary statistics
		print("\nFormatting output")
		split_vcf(args)
		format_vcf(args)
		snps_summary(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory") % utility.max_mem_usage()

	# Optionally remove temporary files
	if args['remove_temp']: remove_tmp(args)

