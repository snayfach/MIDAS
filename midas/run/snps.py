#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, subprocess, shutil, csv
from time import time
from midas import utility
import parse_pileup, Bio.SeqIO

def build_genome_db(args, species):
	""" Build FASTA and BT2 database of representative genomes """
	import Bio.SeqIO
	# fasta database
	outfile = open('/'.join([args['outdir'], 'snps/temp/genomes.fa']), 'w')
	db_stats = {'total_length':0, 'total_seqs':0, 'species':0}
	for sp in species.values():
		db_stats['species'] += 1
		infile = utility.iopen(sp.paths['fna'])
		for r in Bio.SeqIO.parse(infile, 'fasta'):
			outfile.write('>%s\n%s\n' % (r.id, str(r.seq).upper()))
			db_stats['total_length'] += len(r.seq)
			db_stats['total_seqs'] += 1
		infile.close()
	outfile.close()
	# print out database stats
	print("  total genomes: %s" % db_stats['species'])
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
	""" Use Bowtie2 to map reads to representative genomes """
	# Bowtie2
	bam_path = os.path.join(args['outdir'], 'snps/temp/genomes.bam')
	command = '%s --no-unal ' % args['bowtie2']
	command += '-x %s ' % '/'.join([args['outdir'], 'snps/temp/genomes']) # index
	if args['max_reads']: command += '-u %s ' % args['max_reads'] # max num of reads
	if args['trim']: command += '--trim3 %s ' % args['trim'] # trim 3'
	command += '--%s ' % args['speed'] # speed/sensitivity
	command += '--threads %s ' % args['threads'] # threads
	command += '-f ' if args['file_type'] == 'fasta' else '-q ' # input type
	command += '-1 %s -2 %s '  % (args['m1'], args['m2']) if args['m2'] else '-U %s ' % args['m1'] # input reads
	# Pipe to samtools
	command += '| %s view -b - ' % args['samtools'] # convert to bam
	command += '| %s sort -f - %s ' % (args['samtools'], bam_path) # sort bam
	# Run command
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# Check for errors
	utility.check_exit_code(process, command)
	print("  finished aligning")
	print("  checking bamfile integrity")
	utility.check_bamfile(args, bam_path)

def pileup(args):
	""" Filter alignments by % id, use samtools to create pileup, filter low quality bases """
	# Stream bam, filter alignments
	command = 'python %s ' % args['stream_bam']
	command += '%s ' % os.path.join(args['outdir'], 'snps/temp/genomes.bam')
	command += '/dev/stdout '
	command += '%s ' % args['mapid']
	command += '%s ' % args['readq']
	command += '%s ' % args['mapq']
	# Pipe to mpileup
	command += '| %s mpileup '  % args['samtools']
	command += '-d 10000 ' # set max depth
	if not args['baq']: command += '-B ' # BAQ
	if args['adjust_mq']: command += '-C 50 ' # adjust MQ
	if not args['discard']: command += '-A ' # keep discordant read pairs
	command += '-Q %s ' % (args['baseq']) # base quality filtering
	command += '-f %s ' % ('%s/snps/temp/genomes.fa' % args['outdir']) # reference fna file
	command += '- ' #   input bam file
	command += '| gzip > %s ' % ('%s/snps/temp/genomes.mpileup.gz' % args['outdir']) # output file
	# Run command
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)

def write_missing(outfile, ref_id, ref_pos, ref_allele):
	""" Write record for formatted SNP file """
	values = [ref_id, ref_pos, ref_allele, 'NA', '0.0', '0', '0,0,0,0']
	outfile.write('\t'.join(values)+'\n')

def write_present(outfile, pileup):
	""" Write record for formatted SNP file """
	values = [pileup.ref_id, str(pileup.ref_pos), pileup.ref_allele,
	          pileup.alt_allele, str(pileup.ref_freq), str(pileup.depth),
			  pileup.allele_string()]
	outfile.write('\t'.join(values)+'\n')

def format_pileup(args, species, contigs):
	""" Parse mpileups and fill in missing positions """

	# open outfiles
	for sp in species.values():
		sp.out = utility.iopen('/'.join([args['outdir'], 'snps/output/%s.snps.gz' % sp.id]), 'w')
		header = ['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'ref_freq', 'depth', 'count_atcg']
		sp.out.write('\t'.join(header)+'\n')
		sp.contigs = sorted([c.id for c in contigs.values() if c.species_id == sp.id])
		sp.i = 0 # contig index
		sp.j = 0 # position index
		
	# parse pileup
	pileup_path = '%s/snps/temp/genomes.mpileup.gz' % args['outdir']
	for p in parse_pileup.main(pileup_path):
	
		# fetch contig info
		sp = species[contigs[p.ref_id].species_id]
		contig = contigs[sp.contigs[sp.i]]
				
		# contig ids don't match
		#   indicates that one or more upstream contigs have zero coverage
		while p.ref_id != contig.id:
			write_missing(sp.out, ref_id=contig.id, ref_pos=str(sp.j+1), ref_allele=contig.seq[sp.j])
			sp.j += 1
			if sp.j >= contig.length:
				sp.i += 1; sp.j = 0
				contig = contigs[sp.contigs[sp.i]]
			
		# positions don't match
		#   indicates that one or more upstream positions have zero coverage
		while p.ref_pos != sp.j+1:
			write_missing(sp.out, ref_id=contig.id, ref_pos=str(sp.j+1), ref_allele=contig.seq[sp.j])
			sp.j += 1
			if sp.j >= contig.length:
				sp.i += 1; sp.j = 0
				contig = contigs[sp.contigs[sp.i]]
				
		# match
		#   write info from pileup
		write_present(sp.out, pileup=p)
		sp.j += 1
		if sp.j >= contig.length:
			sp.i += 1; sp.j = 0
			contig = contigs[sp.contigs[sp.i]]

	# fill in downstream positions & contigs with zero coverage
	for sp in species.values():
		
		# lefover positions on last contig
		#   indicates that one or more downstream positions have zero coverage
		contig = contigs[sp.contigs[sp.i]]
		while sp.j < contig.length:
			write_missing(sp.out, ref_id=contig.id, ref_pos=str(sp.j+1), ref_allele=contig.seq[sp.j])
			sp.j += 1

		# lefover contigs
		#   indicates that one or more downstream contigs have zero coverage
		sp.i += 1; sp.j = 0
		while sp.i < len(sp.contigs):
			contig = contigs[sp.contigs[sp.i]]
			write_missing(sp.out, ref_id=contig.id, ref_pos=str(sp.j+1), ref_allele=contig.seq[sp.j])
			sp.j += 1
			if sp.j >= contig.length:
				sp.i += 1; sp.j = 0

		# close output files
		sp.out.close()

def snps_summary(args, species):
	""" Get summary of mapping statistics """
	# store stats
	stats = {}
	for species_id in species:
		genome_length, covered_bases, total_depth, maf = [0,0,0,0]
		for r in utility.parse_file('/'.join([args['outdir'], 'snps/output/%s.snps.gz' % species_id])):
			genome_length += 1
			depth = int(r['depth'])
			if depth > 0:
				covered_bases += 1
				total_depth += depth
				ref_freq = float(r['ref_freq'])
				maf += ref_freq if ref_freq <= 0.5 else 1 - ref_freq
		fraction_covered = covered_bases/float(genome_length)
		mean_coverage = total_depth/float(covered_bases) if covered_bases > 0 else 0
		stats[species_id] = {'genome_length':genome_length, 'covered_bases':covered_bases,
							 'fraction_covered':fraction_covered,'mean_coverage':mean_coverage}
	# write stats
	fields = ['genome_length', 'covered_bases', 'fraction_covered', 'mean_coverage']
	outfile = open('/'.join([args['outdir'], 'snps/summary.txt']), 'w')
	outfile.write('\t'.join(['species_id'] + fields)+'\n')
	for species_id in stats:
		record = [species_id] + [str(stats[species_id][field]) for field in fields]
		outfile.write('\t'.join(record)+'\n')
	outfile.close()

def remove_tmp(args):
	""" Remove specified temporary files """
	shutil.rmtree('/'.join([args['outdir'], 'snps/temp']))

class Species:
	""" Base class for species """
	def __init__(self, id):
		self.id = id
		self.paths = {}
		self.contigs = []
		
	def fetch_paths(self, ref_db):
		indir =  '%s/rep_genomes/%s' % (ref_db, self.id)
		for ext in ['', '.gz']:
			for type in ['fna', 'features']:
				path = '%s/genome.%s%s' % (indir, type, ext)
				if os.path.isfile(path):
					self.paths[type] = path

def initialize_species(args):
	species = {}
	splist = '%s/snps/species.txt' % args['outdir']
	if args['build_db']:
		from midas.run.species import select_species
		with open(splist, 'w') as outfile:
			for id in select_species(args):
				species[id] = Species(id)
				outfile.write(id+'\n')
	elif os.path.isfile(splist):
		for line in open(splist):
			id = line.rstrip()
			species[id] = Species(id)
	for sp in species.values():
		sp.fetch_paths(ref_db=args['db'])
	return species

class Contig:
	""" Base class for contig """
	def __init__(self, id):
		self.id = id

def initialize_contigs(species):
	contigs = {}
	for sp in species.values():
		infile = utility.iopen(sp.paths['fna'])
		for rec in Bio.SeqIO.parse(infile, 'fasta'):
			contig = Contig(rec.id)
			contig.seq = str(rec.seq).upper()
			contig.length = len(contig.seq)
			contig.species_id = sp.id
			contigs[contig.id] = contig
		infile.close()
	return contigs

def run_pipeline(args):
	""" Run entire pipeline """
	
	# Initialize reference data
	print("\nReading reference data")
	start = time()
	species = initialize_species(args)
	contigs = initialize_contigs(species)
	print("  %s minutes" % round((time() - start)/60, 2) )
	print("  %s Gb maximum memory" % utility.max_mem_usage())
	
	# Build genome database for selected species
	if args['build_db']:
		print("\nBuilding database of representative genomes")
		args['log'].write("\nBuilding database of representative genomes\n")
		start = time()
		build_genome_db(args, species)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Use bowtie2 to map reads to a representative genome for each species
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

	# Split pileup into files for each species, format, and report summary statistics
		print("\nFormatting output")
		args['log'].write("\nFormatting output\n")
		format_pileup(args, species, contigs)
		snps_summary(args, species)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Optionally remove temporary files
	if args['remove_temp']: remove_tmp(args)

