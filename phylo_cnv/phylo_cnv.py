#!/usr/bin/python

# PhyloCNV - estimating the abundance, gene-content, and phylogeny of microbes from metagenomes
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.1'

# Libraries
# ---------
import sys
import os
import pysam
import subprocess
import resource
from phylo_species import estimate_species_abundance, write_abundance
from gzip import open as gzopen
from Bio.SeqIO import parse as SeqIOparse
from operator import itemgetter
from time import time
from numpy import median
from platform import system
from multiprocessing import Process

# Functions
# ---------

def check_arguments(args):
	""" Check validity of command line arguments """
	
	# Pipeline options
	if not any([args['all'], args['species_profile'],
				args['pangenome_build_db'], args['pangenome_align'], args['pangenome_cov'],
				args['snps_build_db'], args['snps_align'], args['snps_call']]):
		sys.exit('Specify one or more pipeline option(s)')
		
	# Genome cluster selection
	if not any([args['gc_id'], args['gc_topn'], args['gc_cov'], args['gc_rbun']]):
		sys.exit('Specify one or more genome-cluster-selection option(s)')

	# Turn on entire pipeline
	if args['all']:
		args['species_profile'] = True
		args['pangenome_build_db'] = True
		args['pangenome_align'] = True
		args['pangenome_cov'] = True
		args['snps_build_db'] = True
		args['snps_align'] = True
		args['snps_call'] = True

	# Input options
	if not args['m1'] and (args['profile'] or args['align']):
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

def add_paths(args):
	""" Add paths to external files and binaries """
	if system() not in ['Linux', 'Darwin']:
		sys.exit("Operating system '%s' not supported" % system())
	else:
		main_dir = os.path.dirname(os.path.abspath(__file__))
		args['bowtie2-build'] = '/'.join([main_dir, 'bin', system(), 'bowtie2-build'])
		args['bowtie2'] = '/'.join([main_dir, 'bin', system(), 'bowtie2'])
		args['samtools'] = '/'.join([main_dir, 'bin', system(), 'samtools'])
		args['pid_cutoffs'] = '/'.join([main_dir, 'data', 'pid_cutoffs.txt'])
		args['bad_gcs'] = '/'.join([main_dir, 'data', 'bad_cluster_ids.txt'])

def print_copyright():
	# print out copyright information
	print ("-------------------------------------------------------------------------")
	print ("PhyloCNV - estimating the abundance, gene-content, and phylogeny of microbes from metagenomes")
	print ("version %s; github.com/snayfach/PhyloCNV" % __version__)
	print ("Copyright (C) 2015 Stephen Nayfach")
	print ("Freely distributed under the GNU General Public License (GPLv3)")
	print ("-------------------------------------------------------------------------")

def read_phylo_species(inpath):
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

def auto_detect_file_type(inpath):
	""" Detect file type [fasta or fastq] of <p_reads> """
	for line in iopen(inpath):
		if line[0] == '>': return 'fasta'
		elif line[0] == '@': return 'fastq'
		else: sys.exit("Filetype [fasta, fastq] of %s could not be recognized" % inpath)

def iopen(inpath):
	""" Open input file for reading regardless of compression [gzip, bzip] or python version """
	ext = inpath.split('.')[-1]
	# Python2
	if sys.version_info[0] == 2:
		if ext == 'gz': return gzopen(inpath)
		elif ext == 'bz2': return bz2.BZ2File(inpath)
		else: return open(inpath)
	# Python3
	elif sys.version_info[0] == 3:
		if ext == 'gz': return io.TextIOWrapper(gzopen(inpath))
		elif ext == 'bz2': return bz2.BZ2File(inpath)
		else: return open(inpath)

def select_genome_clusters(args):
	""" Select genome clusters to map to """
	cluster_sets = {}
	# read in cluster abundance if necessary
	if any([args['gc_topn'], args['gc_cov'], args['gc_rbun']]):
		inpath = '/'.join([args['out'], 'genome_clusters.abundance'])
		if not os.path.isfile(inpath):
			sys.exit("\nError: genome_clusters.abundance not found in the specified output directory.\nThis file is necessary to specify genome-clusters with --gc_topn, --gc_cov, or --gc_rbun.\nYour options include:\n  1) Include the flag --species_profile when you run phylo_cnv\n  2) Manually specify the genome-clusters you wish to target using --gc_id\n")
		cluster_abundance = read_phylo_species(inpath)
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
		# user specifed a relative-abundance threshold
		if args['gc_topn']:
			cluster_sets['gc_topn'] = set([])
			cluster_abundance = [(i,d['rel_abun']) for i,d in cluster_abundance.items()]
			sorted_abundance = sorted(cluster_abundance, key=itemgetter(1), reverse=True)
			for cluster_id, rel_abun in sorted_abundance[0:args['gc_topn']]:
				cluster_sets['gc_topn'].add(cluster_id)
	# user specified one or more genome-clusters
	if args['gc_id']:
		cluster_sets['gc_rbun'] = set([])
		for cluster_id in args['gc_id']:
			cluster_sets['gc_rbun'].add(cluster_id)
	# intersect sets of genome-clusters
	my_clusters = list(set.intersection(*cluster_sets.values()))
	# check that specified genome-clusters are valid
	for cluster_id in my_clusters:
		if cluster_id not in os.listdir(args['db_dir']):
			sys.exit("\nError: the specified genome_cluster '%s' was not found in the reference database (-D)\n" % cluster_id)
	# remove bad cluster_ids
	for line in open(args['bad_gcs']):
		try: my_clusters.remove(line.rstrip())
		except: pass
	# check that at least one genome-cluster was selected
	if len(my_clusters) == 0:
		sys.exit("\nError: no genome-clusters sastisfied your selection criteria. \n")
	return my_clusters

def pangenome_align(args, tax_mask):
	""" Use Bowtie2 to map reads to all specified genome clusters """
	# Build command
	command = '%s --no-unal ' % args['bowtie2']
	#   index
	command += '-x %s ' % '/'.join([args['out'], 'db', 'pangenomes'])
	#   specify reads
	if args['pangenome_reads']: command += '-u %s ' % args['pangenome_reads']
	#   speed/sensitivity
	command += '--%s ' % args['pangenome_align_speed']
	#   threads
	command += '--threads %s ' % args['threads']
	#   file type
	if args['file_type'] == 'fasta': command += '-f '
	else: command += '-q '
	#   input file
	if (args['m1'] and args['m2']): command += '-1 %s -2 %s ' % (args['m1'], args['m2'])
	else: command += '-U %s' % args['m1']
	#   output unsorted bam
	bampath = '/'.join([args['out'], 'pangenome.bam'])
	command += '| %s view -b - > %s' % (args['samtools'], bampath)
	# Run command
	if args['debug']: print("  running: %s") % command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	#sys.stderr.write(err) # write to stderr: bowtie2 output

def genome_align(args):
	""" Use Bowtie2 to map reads to representative genomes from each genome cluster
	"""
	# Build command
	#	bowtie2
	command = '%s --no-unal ' % args['bowtie2']
	#   index
	command += '-x %s ' % '/'.join([args['out'], 'db', 'genomes'])
	#   specify reads
	if args['snps_reads']: command += '-u %s ' % args['snps_reads']
	#   speed/sensitivity
	command += '--%s ' % args['snps_align_speed']
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
	command += '| %s sort -f - %s ' % (args['samtools'], os.path.join(args['out'], 'genomes.bam'))
	# Run command
	if args['debug']: print("  running: %s") % command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()

def pileup(args):
	""" Use Samtools to create pileup, filter low quality bases, and write results to VCF file """
	# Build command
	#   mpileup
	command = '%s mpileup -uv -A -d 10000 --skip-indels -B ' % args['samtools']
	#   quality filtering
	command += '-q %s -Q %s ' % (args['snps_mapq'], args['snps_baseq'])
	#   reference fna file
	command += '-f %s ' % '/'.join([args['out'], 'db/genomes.fa'])
	#   input bam file
	command += '%s ' % '/'.join([args['out'], 'genomes.bam'])
	#   output vcf file
	command += '> %s ' % '/'.join([args['out'], 'genomes.vcf'])
	# Run command
	if args['debug']: print("  running: %s") % command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
		
def format_vcf(args):
	""" Format vcf output for easy parsing """
	# map scaffold to cluster_id
	ref_to_cluster = {}
	for line in open('/'.join([args['out'], 'db/genomes.map'])):
		ref_id, cluster_id = line.rstrip().split()
		ref_to_cluster[ref_id] = cluster_id
	# open outfiles for each cluster_id
	outdir = '/'.join([args['out'], 'snps'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	outfiles = {}
	genome_clusters = set(ref_to_cluster.values())
	for cluster_id in genome_clusters:
		outfiles[cluster_id] = gzopen('/'.join([outdir, '%s.snps.gz' % cluster_id]), 'w')
		outfiles[cluster_id].write('\t'.join(['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
											  'count_alleles', 'count_ref', 'count_alt', 'depth', 'ref_freq'])+'\n')
	# parse vcf into snp files for each cluster_id
	inpath = '/'.join([args['out'], 'genomes.vcf'])
	for r in parse_vcf(inpath):
		rec = [r['ref_id'], r['ref_pos'], r['ref_allele'], r['alt_allele'], r['cons_allele'],
			   r['count_alleles'], r['count_ref'], r['count_alt'], r['depth'], r['ref_freq']]
		outfile = outfiles[ref_to_cluster[r['ref_id']]]
		outfile.write('\t'.join([str(x) for x in rec])+'\n')

def parse_vcf(inpath):
	""" Yields formatted records from VCF output """
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

def compute_perc_id(aln):
	""" Compute percent identity for read """
	length = aln.query_length
	edit = dict(aln.tags)['NM']
	return 100 * (length - edit)/float(length)

def compute_aln_cov(aln):
	""" Compute percent identity for paired-end read """
	aln_cov = len(aln.query_alignment_sequence)/float(aln.query_length)
	return aln_cov

def count_mapped_bp(args):
	""" Count number of bp mapped to each centroid across pangenomes """
	bam_path = '/'.join([args['out'], 'pangenome.bam'])
	aln_file = pysam.AlignmentFile(bam_path, "rb")
	ref_to_length = dict([(i,j) for i,j in zip(aln_file.references, aln_file.lengths)])
	ref_to_cov = dict([(i,0.0) for i in aln_file.references])
	for aln in aln_file.fetch(until_eof = True):
		query = aln.query_name
		if compute_perc_id(aln) < args['pangenome_map_pid']:
			continue
		elif compute_aln_cov(aln) < args['pangenome_aln_cov']:
			continue
		else:
			ref_id = aln_file.getrname(aln.reference_id)
			cov = len(aln.query_alignment_sequence)/float(ref_to_length[ref_id])
			ref_to_cov[ref_id] += cov
	return ref_to_cov

def read_centroid_map(args, genome_clusters):
	""" Map 99% ID pangenome centroids to a lower level (90, 92.5, 95, 97.5)"""
	centroid_map = {}
	for cluster_id in genome_clusters:
		inpath = '/'.join([args['db_dir'], cluster_id, 'gene_family_map.txt.gz'])
		infile = gzopen(inpath)
		next(infile)
		for line in infile:
			x = line.rstrip().split()
			y = {'90':x[0], '92.5':x[1], '95':x[2], '97.5':x[3], '99':x[4]}
			centroid_map[y['99']] = y[args['pangenome_pid']]
	return centroid_map

def aggregate_coverage(args, ref_to_cov, genome_clusters):
	""" Aggregate gene coverage at given percent identity clustering """
	centroid_map = read_centroid_map(args, genome_clusters)

	centroid_to_cov = dict([(x, 0.0) for x in set(centroid_map.values())])
	for ref_id, cov in ref_to_cov.items():
		ref_id2 = centroid_map[ref_id]
		centroid_to_cov[ref_id2] += cov
	return centroid_to_cov

def compute_phyeco_cov(args, genome_clusters):
	""" Count number of bp mapped to each PhyEco marker gene """
	# read in cutoffs
	phyeco_to_cutoff = {}
	for line in open(args['pid_cutoffs']):
		phyeco_id, pid = line.rstrip().split()
		phyeco_to_cutoff[phyeco_id] = float(pid)
	# read in map of gene to phyeco
	ref_to_cluster = {}
	ref_to_phyeco = {}
	for cluster_id in genome_clusters:
		inpath = '/'.join([args['db_dir'], cluster_id, 'universal_genes.txt.gz'])
		infile = gzopen(inpath)
		next(infile)
		for line in infile:
			gene_id, phyeco_id = line.rstrip().split()
			ref_to_phyeco[gene_id] = phyeco_id
			ref_to_cluster[gene_id] = cluster_id
	# init phyeco coverage
	cluster_to_phyeco_to_cov = {}
	for cluster_id in genome_clusters:
		cluster_to_phyeco_to_cov[cluster_id] = {}
		for phyeco_id in phyeco_to_cutoff:
			cluster_to_phyeco_to_cov[cluster_id][phyeco_id] = 0.0
	# open input file
	bam_path = '/'.join([args['out'], 'pangenome.bam'])
	aln_file = pysam.AlignmentFile(bam_path, "rb")
	# read in map of gene lengths
	ref_to_length = dict([(i,j) for i,j in zip(aln_file.references, aln_file.lengths)])
	# count bp mapped to phyeco genes
	for aln in aln_file.fetch(until_eof = True):
		ref_id = aln_file.getrname(aln.reference_id)
		phyeco_id = ref_to_phyeco[ref_id] if ref_id in ref_to_phyeco else None
		if not phyeco_id:
			continue
		elif phyeco_id not in phyeco_to_cutoff:
			continue
		elif compute_perc_id(aln) < phyeco_to_cutoff[phyeco_id]:
			continue
		elif compute_aln_cov(aln) < args['pangenome_aln_cov']:
			continue
		else:
			cluster_id = ref_to_cluster[ref_id]
			cov = len(aln.query_alignment_sequence)/float(ref_to_length[ref_id])
			cluster_to_phyeco_to_cov[cluster_id][phyeco_id] += cov
	# compute median phyeco cov
	cluster_to_norm = {}
	for cluster_id in cluster_to_phyeco_to_cov:
		covs = cluster_to_phyeco_to_cov[cluster_id].values()
		cluster_to_norm[cluster_id] = median(covs)
	return cluster_to_norm

def compute_pangenome_coverage(args):
	""" Compute coverage of pangenome for cluster_id and write results to disk """
	# map ref_id to cluster_id
	ref_to_cluster = {}
	for line in open('/'.join([args['out'], 'db/pangenome.map'])):
		ref_id, cluster_id = line.rstrip().split()
		ref_to_cluster[ref_id] = cluster_id
	# open outfiles for each cluster_id
	outdir = '/'.join([args['out'], 'coverage'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	outfiles = {}
	genome_clusters = set(ref_to_cluster.values())
	for cluster_id in genome_clusters:
		outfiles[cluster_id] = gzopen('/'.join([outdir, '%s.cov.gz' % cluster_id]), 'w')
		outfiles[cluster_id].write('\t'.join(['ref_id', 'raw_coverage', 'normalized_coverage'])+'\n')
	# parse bam into cov files for each cluster_id
	cluster_to_norm = compute_phyeco_cov(args, genome_clusters)
	ref_to_cov = count_mapped_bp(args)
	# aggregate coverages
	if args['pangenome_pid'] != '99':
		ref_to_cov = aggregate_coverage(args, ref_to_cov, genome_clusters)
	# write to output files
	for ref_id, cov in ref_to_cov.items():
		outfile = outfiles[ref_to_cluster[ref_id]]
		normcov = cov/cluster_to_norm[cluster_id] if cluster_to_norm[cluster_id] > 0 else 0
		outfile.write('\t'.join([str(x) for x in [ref_id, cov, normcov]])+'\n')

def max_mem_usage():
	""" Return max mem usage (Gb) of self and child processes """
	max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
	return round((max_mem_self + max_mem_child)/float(1e6), 2)

def build_pangenome_db(args, genome_clusters):
	""" Build FASTA and BT2 database from pangene cluster centroids """
	# fasta database
	outdir = '/'.join([args['out'], 'db'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	pangenome_fasta = open('/'.join([args['out'], 'db/pangenomes.fa']), 'w')
	pangenome_map = open('/'.join([args['out'], 'db/pangenome.map']), 'w')
	db_stats = {'total_length':0, 'total_seqs':0, 'genome_clusters':0}
	for cluster_id in genome_clusters:
		db_stats['genome_clusters'] += 1
		inpath = '/'.join([args['db_dir'], cluster_id, 'pangenome.fa.gz'])
		infile = gzopen(inpath)
		for r in SeqIOparse(infile, 'fasta'):
			genome_id = '.'.join(r.id.split('.')[0:2])
			if not args['tax_mask'] or genome_id not in args['tax_mask']:
				pangenome_fasta.write('>%s\n%s\n' % (r.id, str(r.seq)))
				pangenome_map.write('%s\t%s\n' % (r.id, cluster_id))
				db_stats['total_length'] += len(r.seq)
				db_stats['total_seqs'] += 1
	# print out database stats
	if args['verbose']:
		print("  total genome-clusters: %s" % db_stats['genome_clusters'])
		print("  total genes: %s" % db_stats['total_seqs'])
		print("  total base-pairs: %s" % db_stats['total_length'])
	# bowtie2 database
	inpath = '/'.join([args['out'], 'db/pangenomes.fa'])
	outpath = '/'.join([args['out'], 'db/pangenomes'])
	command = ' '.join([args['bowtie2-build'], inpath, outpath])
	if args['debug']: print("  running: %s" % command)
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()

def fetch_centroid(args, cluster_id):
	""" Get the genome_id corresponding to cluster centroid """
	inpath = '/'.join([args['db_dir'], cluster_id, 'genomes.txt.gz'])
	infile = gzopen(inpath)
	for line in infile:
		if line.split()[2] == 'Y':
			return line.split()[1]

def build_genome_db(args, genome_clusters):
	""" Build FASTA and BT2 database from genome cluster centroids """
	# fasta database
	outdir = '/'.join([args['out'], 'db'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	genomes_fasta = open('/'.join([args['out'], 'db', 'genomes.fa']), 'w')
	genomes_map = open('/'.join([args['out'], 'db', 'genomes.map']), 'w')
	db_stats = {'total_length':0, 'total_seqs':0, 'genome_clusters':0}
	for cluster_id in genome_clusters:
		if args['tax_mask'] and fetch_centroid(args, cluster_id) in args['tax_mask']:
			continue
		db_stats['genome_clusters'] += 1
		inpath = '/'.join([args['db_dir'], cluster_id, 'representative.fna.gz'])
		infile = gzopen(inpath)
		for line in infile:
			genomes_fasta.write(line)
			db_stats['total_length'] += len(line.rstrip())
			if line[0] == '>':
				sid = line.rstrip().lstrip('>').split()[0]
				genomes_map.write(sid+'\t'+cluster_id+'\n')
				db_stats['total_seqs'] += 1
	# print out database stats
	if args['verbose']:
		print("  total genomes: %s" % db_stats['genome_clusters'])
		print("  total contigs: %s" % db_stats['total_seqs'])
		print("  total base-pairs: %s" % db_stats['total_length'])
	# bowtie2 database
	inpath = '/'.join([args['out'], 'db', 'genomes.fa'])
	outpath = '/'.join([args['out'], 'db', 'genomes'])
	command = ' '.join([args['bowtie2-build'], inpath, outpath])
	if args['debug']: print("  running: %s" % command)
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()

def run_pipeline(args):
	""" Run entire pipeline """
	
	check_arguments(args) # Check validity of command line arguments
	add_paths(args) # Add paths to external files and binaries
	if args['verbose']: print_copyright()
	args['file_type'] = auto_detect_file_type(args['m1'])

	# 1a. Estimate the abundance of genome-clusters
	if args['species_profile']:
		start = time()
		if args['verbose']: print("\nEstimating the abundance of genome-clusters")
		cluster_abundance = estimate_species_abundance({'inpath':args['m1'],'nreads':args['reads_ms'],'threads':args['threads']})
		write_abundance('%s/genome_clusters.abundance' % args['out'], cluster_abundance)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 2a. Build pangenome database for selected GCs
	if args['pangenome_build_db']:
		if args['verbose']: print("\nBuilding pangenome database")
		start = time()
		genome_clusters = select_genome_clusters(args)
		build_pangenome_db(args, genome_clusters)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 2a. Use bowtie2 to align reads to pangenome database
	if args['pangenome_align']:
		start = time()
		if args['verbose']: print("\nAligning reads to pangenomes")
		pangenome_align(args, args['tax_mask'])
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 2b. Compute pangenome coverage for each genome-cluster
	if args['pangenome_cov']:
		start = time()
		if args['verbose']: print("\nComputing coverage of pangenomes")
		compute_pangenome_coverage(args)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 3a. Build genome database for selected GCs
	if args['snps_build_db']:
		if args['verbose']: print("\nBuilding database of representative genomes")
		start = time()
		genome_clusters = select_genome_clusters(args)
		build_genome_db(args, genome_clusters)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 3b. Use bowtie2 to map reads to a representative genome for each genome-cluster
	if args['snps_align']:
		if args['verbose']: print("\nMapping reads to representative genomes")
		start = time()
		genome_align(args)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 3c. Use mpileup to identify SNPs
	if args['snps_call']:
		start = time()
		if args['verbose']: print("\nRunning mpileup")
		pileup(args)
		format_vcf(args)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

