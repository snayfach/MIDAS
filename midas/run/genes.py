#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, subprocess, gzip, csv, Bio.SeqIO, numpy as np
from collections import defaultdict
from time import time
from midas import utility

class Species:
	""" Base class for species """
	def __init__(self, id):
		self.id = id
		self.paths = {}
		self.genes = []
		self.pangenome_size = 0
		self.aligned_reads = 0
		self.mapped_reads = 0
		self.markers = defaultdict(float)
		
	def init_ref_db(self, ref_db):
		""" Set paths to input files """
		self.dir = '%s/pan_genomes/%s' % (ref_db, self.id)
		for ext in ['', '.gz']:
			for file in ['centroids.ffn', 'cluster_info.txt', 'gene_info.txt']:
				inpath = '%s/%s%s' % (self.dir, file, ext)
				if os.path.isfile(inpath):
					self.paths[file] = inpath

def initialize_species(args):
	""" Initialize Species objects """
	species = {}
	splist = '%s/genes/species.txt' % args['outdir']
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
		sp.init_ref_db(args['db'])
	return species

class Gene:
	""" Base class for gene """
	def __init__(self, id):
		self.id = id
		self.aligned_reads = 0
		self.mapped_reads = 0
		self.depth = 0.0
		self.length = 0
		self.copies = 0.0
		self.marker_id = None

def initialize_genes(args, species):
	""" Initialize Gene objects """
	genes = {}
	# fetch gene_id, species_id, gene length
	for sp in species.values():
		path = sp.paths['centroids.ffn']
		file = utility.iopen(path)
		for seq in Bio.SeqIO.parse(file, 'fasta'):
			genes[seq.id] = Gene(seq.id)
			genes[seq.id].species_id = sp.id
			genes[seq.id].length = len(seq.seq)
			sp.pangenome_size += 1
		file.close()
	# fetch marker_id
	path = '%s/marker_genes/phyeco.map' % args['db']
	file = utility.iopen(path)
	reader = csv.DictReader(file, delimiter='\t')
	for r in reader:
		if r['gene_id'] in genes:
			genes[r['gene_id']].marker_id=r['marker_id']
	file.close()
	return genes

def build_pangenome_db(args, species):
	""" Build FASTA and BT2 database from pangene species centroids """
	import Bio.SeqIO
	# fasta database
	outdir = '/'.join([args['outdir'], 'genes/temp'])
	pangenome_fasta = open('/'.join([outdir, 'pangenomes.fa']), 'w')
	pangenome_map = open('/'.join([outdir, 'pangenomes.map']), 'w')
	db_stats = {'total_length':0, 'total_seqs':0, 'species':0}
	for sp in species.values():
		db_stats['species'] += 1
		infile = utility.iopen(sp.paths['centroids.ffn'])
		for r in Bio.SeqIO.parse(infile, 'fasta'):
			pangenome_fasta.write('>%s\n%s\n' % (r.id, str(r.seq).upper()))
			pangenome_map.write('%s\t%s\n' % (r.id, sp.id))
			db_stats['total_length'] += len(r.seq)
			db_stats['total_seqs'] += 1
		infile.close()
	pangenome_fasta.close()
	pangenome_map.close()
	# print out database stats
	print("  total species: %s" % db_stats['species'])
	print("  total genes: %s" % db_stats['total_seqs'])
	print("  total base-pairs: %s" % db_stats['total_length'])
	# bowtie2 database
	command = '%s ' % args['bowtie2-build']
	command += '--threads %s ' % args['threads']
	command += '%s/genes/temp/pangenomes.fa ' % args['outdir']
	command += '%s/genes/temp/pangenomes ' % args['outdir']
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)

def pangenome_align(args):
	""" Use Bowtie2 to map reads to all specified genome species """
	# Bowtie2
	command = '%s --no-unal ' % args['bowtie2']
	command += '-x %s ' % '/'.join([args['outdir'], 'genes/temp/pangenomes']) # index
	if args['max_reads']: command += '-u %s ' % args['max_reads'] # max num of reads
	if args['trim']: command += '--trim3 %s ' % args['trim'] # trim 3'
	command += '--%s-local ' % args['speed'] #   speed/sensitivity
	command += '--threads %s ' % args['threads'] #   threads
	command += '-f ' if args['file_type'] == 'fasta' else '-q ' # input type
	if args['m2']: # -1 and -2 contain paired reads
		command += '-1 %s -2 %s ' % (args['m1'], args['m2']) 
	elif args['interleaved']: # -1 contains paired reads
		command += '--interleaved %s ' % args['m1'] 
	else: # -1 contains unpaired reads
		command += '-U %s ' % args['m1'] 
	# Output unsorted bam
	bampath = '/'.join([args['outdir'], 'genes/temp/pangenomes.bam'])
	command += '| %s view ' % args['samtools']
	command += '--threads %s ' % args['threads']
	command += '-b - > %s' % bampath
	# Run command
	args['log'].write('command: '+command+'\n')
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# Check for errors
	utility.check_exit_code(process, command)
	print("  finished aligning")
	print("  checking bamfile integrity")
	utility.check_bamfile(args, bampath)

def pangenome_coverage(args, species, genes):
	""" Compute coverage of pangenome for species_id and write results to disk """
	count_mapped_bp(args, species, genes)
	normalize(args, species, genes)
	write_results(args, species, genes)

def keep_read(aln, min_pid, min_readq, min_mapq, min_aln_cov):
	align_len = len(aln.query_alignment_sequence)
	query_len = aln.query_length
	# min pid
	if 100*(align_len-dict(aln.tags)['NM'])/float(align_len) < min_pid:
		return False
	# min read quality
	elif np.mean(aln.query_qualities) < min_readq:
		return False
	# min map quality
	elif aln.mapping_quality < min_mapq:
		return False
	# min aln cov
	elif align_len/float(query_len)  < min_aln_cov:
		return False
	else:
		return True
		
def count_mapped_bp(args, species, genes):
	""" Count number of bp mapped to each gene across pangenomes """
	import pysam
	bam_path = '/'.join([args['outdir'], 'genes/temp/pangenomes.bam'])
	bamfile = pysam.AlignmentFile(bam_path, "rb")

	# loop over alignments, sum values per gene
	for index, aln in enumerate(bamfile.fetch(until_eof = True)):
			
		gene = genes[bamfile.getrname(aln.reference_id)]
		species[gene.species_id].aligned_reads += 1
		gene.aligned_reads += 1
		
		if not keep_read(aln, args['mapid'], args['readq'], args['mapq'], args['aln_cov']):
			continue
		else:
			species[gene.species_id].mapped_reads += 1			
			gene.mapped_reads += 1
			gene.depth += len(aln.query_alignment_sequence)/float(gene.length)
	
	print("  total aligned reads: %s" % sum([sp.aligned_reads for sp in species.values()]))
	print("  total mapped reads: %s" % sum([sp.mapped_reads for sp in species.values()]))
	
	# loop over genes, sum values per species
	for gene in genes.values():
		species[gene.species_id].genes.append(gene.depth)
	
	# loop over species, compute summaries
	for sp in species.values():
		non_zero_genes = [_ for _ in sp.genes if _ > 0]
		sp.covered_genes = len(non_zero_genes)
		sp.mean_coverage = np.mean(non_zero_genes) if len(non_zero_genes) > 0 else 0
		sp.fraction_covered = sp.covered_genes/float(sp.pangenome_size)

def normalize(args, species, genes):
	""" Count number of bp mapped to each marker gene """
	# compute marker depth
	for gene in genes.values():
		if gene.marker_id is not None:
			species[gene.species_id].markers[gene.marker_id] += gene.depth
	# compute median marker depth
	for sp in species.values():
		sp.marker_coverage = np.median(list(sp.markers.values()))
	# normalize genes by median marker depth
	for gene in genes.values():
		sp = species[gene.species_id]
		if sp.marker_coverage > 0:
			gene.copies = gene.depth/sp.marker_coverage

def write_results(args, species, genes):
	""" Write results to disk """
	# open outfiles for each species_id
	header = ['gene_id', 'count_reads', 'coverage', 'copy_number']
	for sp in species.values():
		path = '/'.join([args['outdir'], 'genes/output/%s.genes.gz' % sp.id])
		sp.out = utility.iopen(path, 'w')
		sp.out.write('\t'.join(header)+'\n')
	# write to output files
	for gene_id in sorted(genes):
		gene = genes[gene_id]
		sp = species[gene.species_id]
		values = [gene.id, gene.mapped_reads, gene.depth, gene.copies]
		sp.out.write('\t'.join([str(_) for _ in values])+'\n')
	# close output files
	for sp in species.values():
		sp.out.close()
	# summary stats
	path = '/'.join([args['outdir'], 'genes/summary.txt'])
	file = open(path, 'w')
	header = ['species_id', 'pangenome_size', 'covered_genes', 'fraction_covered', 'mean_coverage', 'marker_coverage', 'aligned_reads', 'mapped_reads']
	file.write('\t'.join(header)+'\n')
	for sp in species.values():
		values = [sp.id, sp.pangenome_size, sp.covered_genes, sp.fraction_covered, sp.mean_coverage, sp.marker_coverage, sp.aligned_reads, sp.mapped_reads]
		file.write('\t'.join([str(_) for _ in values])+'\n')
	file.close()

def remove_tmp(args):
	""" Remove specified temporary files """
	import shutil
	shutil.rmtree('/'.join([args['outdir'], 'genes/temp']))

def run_pipeline(args):
	""" Run entire pipeline """
	
	# Initialize reference data
	print("\nReading reference data")
	start = time()
	species = initialize_species(args)
	genes = initialize_genes(args, species)
	print("  %s minutes" % round((time() - start)/60, 2) )
	print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Build pangenome database for selected species
	if args['build_db']:
		print("\nBuilding pangenome database")
		args['log'].write("\nBuilding pangenome database\n")
		start = time()
		build_pangenome_db(args, species)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Use bowtie2 to align reads to pangenome database
	if args['align']:
		start = time()
		print("\nAligning reads to pangenomes")
		args['log'].write("\nAligning reads to pangenomes\n")
		pangenome_align(args)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Compute pangenome coverage for each species
	if args['cov']:
		start = time()
		print("\nComputing coverage of pangenomes")
		args['log'].write("\nComputing coverage of pangenomes\n")
		pangenome_coverage(args, species, genes)
		print("  %s minutes" % round((time() - start)/60, 2) )
		print("  %s Gb maximum memory" % utility.max_mem_usage())

	# Optionally remove temporary files
	if args['remove_temp']: remove_tmp(args)
		


