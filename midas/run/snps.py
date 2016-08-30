#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, subprocess, gzip
from time import time
from midas import utility

def build_genome_db(args, genome_clusters):
	""" Build FASTA and BT2 database from genome cluster centroids """
	# fasta database
	genomes_fasta = open('/'.join([args['outdir'], 'snps/temp/genomes.fa']), 'w')
	genomes_map = open('/'.join([args['outdir'], 'snps/temp/genomes.map']), 'w')
	db_stats = {'total_length':0, 'total_seqs':0, 'genome_clusters':0}
	for species_id in genome_clusters:
		if args['tax_mask'] and fetch_centroid(args, species_id) in args['tax_mask']:
			continue
		db_stats['genome_clusters'] += 1
		inpath = '/'.join([args['db'], 'genome_clusters', species_id, 'genome.fna.gz'])
		infile = utility.iopen(inpath)
		for line in infile:
			genomes_fasta.write(line)
			db_stats['total_length'] += len(line.rstrip())
			if line[0] == '>':
				sid = line.rstrip().lstrip('>').split()[0]
				genomes_map.write(sid+'\t'+species_id+'\n')
				db_stats['total_seqs'] += 1
	# print out database stats
	print("  total genomes: %s" % db_stats['genome_clusters'])
	print("  total contigs: %s" % db_stats['total_seqs'])
	print("  total base-pairs: %s" % db_stats['total_length'])
	# bowtie2 database
	inpath = '/'.join([args['outdir'], 'snps/temp/genomes.fa'])
	outpath = '/'.join([args['outdir'], 'snps/temp/genomes'])
	command = ' '.join([args['bowtie2-build'], inpath, outpath])
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)

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
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# Check for errors
	print("   finished aligning")
	print("   checking bamfile integrity")
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
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)

def split_vcf(args):
	""" Format vcf output for easy parsing """
	inpath = os.path.join(args['outdir'], 'snps/temp/genomes.map')
	ref_to_species = utility.read_ref_to_cluster(inpath)
	# open outfiles for each species_id
	outdir = '/'.join([args['outdir'], 'snps/temp/vcf'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	outfiles = {}
	for species_id in set(ref_to_species.values()):
		outfiles[species_id] = open('/'.join([outdir, '%s.vcf' % species_id]), 'w')
	# parse vcf into temorary vcf files for each species_id
	for line in open('/'.join([args['outdir'], 'snps/temp/genomes.vcf'])):
		if line[0] == '#': continue
		species_id = ref_to_species[line.split()[0]]
		outfiles[species_id].write(line)
	# close outfiles
	for file in outfiles.values():
		file.close()

def read_ref_bases(args, species_id):
	""" Read in reference genome by position """
	import Bio.SeqIO
	ref = []
	centroid_path = '/'.join([args['db'], 'genome_clusters', species_id, 'genome.fna.gz'])
	infile = utility.iopen(centroid_path)
	for rec in Bio.SeqIO.parse(infile, 'fasta'):
		for pos in range(1, len(rec.seq)+1):
			ref.append([rec.id, pos, rec.seq[pos-1].upper()])
	return sorted(ref)

def write_snp_record(outfile, snp=None, ref=None, header=False):
	""" Write record for formatted SNP file """
	fields = ['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
	          'count_ref', 'count_alt', 'depth', 'ref_freq']
	if header: # just write header
		outfile.write('\t'.join(fields)+'\n')
	elif ref: # missing snp
		snp = {'ref_id': ref[0], 'ref_pos': str(ref[1]), 'ref_allele': ref[2], 'alt_allele': 'NA',
			   'cons_allele': 'NA', 'depth': 0, 'count_ref': 0, 'count_alt': 0, 'ref_freq': 'NA'}
		record = [str(snp[field]) for field in fields]
		outfile.write('\t'.join(record)+'\n')
	else: # present snp
		record = [str(snp[field]) for field in fields]
		outfile.write('\t'.join(record)+'\n')

def format_vcf(args):
	""" Format vcf files to snp files and fill in missing positions """
	inpath = os.path.join(args['outdir'], 'snps/temp/genomes.map')
	ref_to_species = utility.read_ref_to_cluster(inpath)
	for species_id in set(ref_to_species.values()):
		# open outfile
		outpath = '/'.join([args['outdir'], 'snps/output/%s.snps.gz' % species_id])
		outfile = utility.iopen(outpath, 'w')
		write_snp_record(outfile, header=True)
		# read sorted reference
		ref = read_ref_bases(args, species_id)
		ref_index = 0
		ref_length = len(ref)
		# write formatted records
		vcf_path = '/'.join([args['outdir'], 'snps/temp/vcf/%s.vcf' % species_id])
		for snp in parse_vcf(vcf_path): # loop over formatted records from vcf
			while [snp['ref_id'], snp['ref_pos']] != ref[ref_index][0:2]: # fill in missing snp positions
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
			   'ref_pos':int(r[1]),
			   'site_id':r[0]+'|'+r[1],
			   'ref_allele':r[3],
			   'alt_allele':alt_alleles[0] if len(alt_alleles) > 0 else 'NA',
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
	ref_to_species = utility.read_ref_to_cluster(inpath)
	for species_id in set(ref_to_species.values()):
		genome_length, covered_bases, total_depth, identity, maf = [0,0,0,0,0]
		for r in utility.parse_file('/'.join([args['outdir'], 'snps/output/%s.snps.gz' % species_id])):
			genome_length += 1
			depth = int(r['depth'])
			if depth > 0:
				covered_bases += 1
				total_depth += depth
				if r['ref_allele'] == r['cons_allele']:
					identity += 1
				ref_freq = float(r['ref_freq'])
				maf += ref_freq if ref_freq <= 0.5 else 1 - ref_freq
		stats[species_id] = {'genome_length':genome_length,
							 'covered_bases': covered_bases,
							 'fraction_covered':covered_bases/float(genome_length),
							 'mean_coverage':total_depth/float(covered_bases) if covered_bases > 0 else 0
							 }
	# write stats
	fields = ['genome_length', 'covered_bases', 'fraction_covered', 'mean_coverage']
	outfile = open('/'.join([args['outdir'], 'snps/summary.txt']), 'w')
	outfile.write('\t'.join(['species_id'] + fields)+'\n')
	for species_id in stats:
		record = [species_id] + [str(stats[species_id][field]) for field in fields]
		outfile.write('\t'.join(record)+'\n')

def fetch_centroid(args, species_id):
	""" Get the genome_id corresponding to cluster centroid """
	inpath = '/'.join([args['db'], 'genome_clusters', species_id, 'genomes.txt.gz'])
	infile = utility.iopen(inpath)
	for line in infile:
		if line.split()[2] == 'Y':
			return line.split()[1]

def remove_tmp(args):
	""" Remove specified temporary files """
	import shutil
	shutil.rmtree('/'.join([args['outdir'], 'snps/temp']))

def run_pipeline(args):
	""" Run entire pipeline """
	
	# Build genome database for selected GCs
	if args['build_db']:
		from midas.run import species
		print("\nBuilding database of representative genomes")
		args['log'].write("\nBuilding database of representative genomes\n")
		start = time()
		genome_clusters = species.select_genome_clusters(args)
		build_genome_db(args, genome_clusters)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Use bowtie2 to map reads to a representative genome for each genome-cluster
	if args['align']:
		args['file_type'] = utility.auto_detect_file_type(args['m1'])
		print("\nMapping reads to representative genomes")
		args['log'].write("\nMapping reads to representative genomes\n")
		start = time()
		genome_align(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Use mpileup to identify SNPs
	if args['call']:
		start = time()
		print("\nRunning mpileup")
		args['log'].write("\nRunning mpileup\n")
		pileup(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Split vcf into files for each GC, format, and report summary statistics
		print("\nFormatting output")
		args['log'].write("\nFormatting output\n")
		split_vcf(args)
		format_vcf(args)
		snps_summary(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Optionally remove temporary files
	if args['remove_temp']: remove_tmp(args)

