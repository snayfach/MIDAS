#!/usr/bin/python

# PhyloCNV - estimating the abundance, gene-content, and phylogeny of microbes from metagenomes
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.2'

# Libraries
# ---------
import sys
import os
import pysam
import subprocess
import resource
import gzip
import Bio.SeqIO
import operator
from time import time
from numpy import median
from platform import system
from multiprocessing import Process
from phylo_species import estimate_species_abundance, write_abundance


# Functions
# ---------

def check_arguments(args):
	""" Check validity of command line arguments """
	
	# Pipeline options
	if not any([args['all'], args['species_profile'],
				args['pangenome_build_db'], args['pangenome_align'], args['pangenome_cov'],
				args['snps_build_db'], args['snps_align'], args['snps_call']]):
		sys.exit('Specify one or more pipeline option(s)')
		
	# Turn on entire pipeline
	if args['all']:
		args['species_profile'] = True
		args['pangenome_build_db'] = True
		args['pangenome_align'] = True
		args['pangenome_cov'] = True
		args['snps_build_db'] = True
		args['snps_align'] = True
		args['snps_call'] = True

	# Genome cluster selection
	build_db = (args['pangenome_build_db'] or args['snps_build_db'])
	select_gc = any([args['gc_id'], args['gc_topn'], args['gc_cov'], args['gc_rbun']])
	if build_db and not select_gc:
		sys.exit('To build a reference database, you must specify genome-clusters: --gc_id, --gc_topn, --gc_cov, and/or --gc_rbun')

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
		args['filter_bam'] = '/'.join([main_dir, 'bin', 'filter_bam.py'])

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
		if ext == 'gz': return gzip.open(inpath)
		elif ext == 'bz2': return bz2.BZ2File(inpath)
		else: return open(inpath)
	# Python3
	elif sys.version_info[0] == 3:
		if ext == 'gz': return io.TextIOWrapper(gzip.open(inpath))
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
			sorted_abundance = sorted(cluster_abundance, key=operator.itemgetter(1), reverse=True)
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
	""" Filter alignments by % id, use samtools to create pileup, filter low quality bases, and write results to VCF file """
	# Build command
	#   percent id filtering
	command  = 'python %s %s %s %s | ' % (args['filter_bam'], '%s/genomes.bam' % args['out'], '/dev/stdout', args['snps_pid'])
	#   mpileup
	command += '%s mpileup -uv -A -d 10000 --skip-indels -B ' % args['samtools']
	#   quality filtering
	command += '-q %s -Q %s ' % (args['snps_mapq'], args['snps_baseq'])
	#   reference fna file
	command += '-f %s ' % ('%s/db/genomes.fa' % args['out'])
	#   input bam file
	command += '- '
	#   output vcf file
	command += '> %s ' % ('%s/genomes.vcf' % args['out'])
	# Run command
	if args['debug']: print("  running: %s") % command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()

def read_ref_to_cluster(args):
	""" Read in map of scaffold id to genome-cluster id """
	ref_to_cluster = {}
	for line in open('/'.join([args['out'], 'db/genomes.map'])):
		ref_id, cluster_id = line.rstrip().split()
		ref_to_cluster[ref_id] = cluster_id
	return ref_to_cluster

def split_vcf(args):
	""" Format vcf output for easy parsing """
	ref_to_cluster = read_ref_to_cluster(args)
	# open outfiles for each cluster_id
	outdir = '/'.join([args['out'], 'vcf'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	outfiles = {}
	for cluster_id in set(ref_to_cluster.values()):
		outfiles[cluster_id] = open('/'.join([outdir, '%s.vcf' % cluster_id]), 'w')
	# parse vcf into temorary vcf files for each cluster_id
	for line in open('/'.join([args['out'], 'genomes.vcf'])):
		if line[0] == '#': continue
		cluster_id = ref_to_cluster[line.split()[0]]
		outfiles[cluster_id].write(line)
	# close outfiles
	for file in outfiles.values():
		file.close()

def read_ref_bases(args, cluster_id):
	""" Read in reference genome by position """
	ref = []
	centroid_path = '/'.join([args['db_dir'],cluster_id,'representative.fna.gz'])
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
	outdir = '/'.join([args['out'], 'snps'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	for cluster_id in set(read_ref_to_cluster(args).values()):
		# open outfile
		outfile = gzip.open('/'.join([outdir, '%s.snps.gz' % cluster_id]), 'w')
		write_snp_header(outfile)
		# read sorted reference
		ref = read_ref_bases(args, cluster_id)
		ref_index = 0
		ref_length = len(ref)
		# write formatted records
		vcf_path = '/'.join([args['out'], 'vcf', '%s.vcf' % cluster_id])
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

def parse_snps(inpath):
	""" Yields formatted records from SNPs output """
	infile = gzip.open(inpath)
	next(infile)
	fields = ['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
			  'count_alleles', 'count_ref', 'count_alt', 'depth', 'ref_freq']
	for line in infile:
		yield dict([(i,j) for i,j in zip(fields, line.rstrip().split())])

def snps_summary(args):
	""" Get summary of mapping statistics """
	# store stats
	stats = {}
	for cluster_id in set(read_ref_to_cluster(args).values()):
		genome_length, covered_bases, total_depth, identity, maf = [0,0,0,0,0]
		for r in parse_snps('/'.join([args['out'], 'snps/%s.snps.gz' % cluster_id])):
			genome_length += 1
			depth = int(r['depth'])
			if depth > 0:
				covered_bases += 1
				total_depth += depth
				if r['ref_allele'] == r['cons_allele']:
					identity += 1
				ref_freq = float(r['ref_freq'])
				maf += ref_freq if ref_freq <= 0.5 else 1 - ref_freq
		stats[cluster_id] = {'genome_length':genome_length, 'fraction_covered':covered_bases/float(genome_length),
							 'average_depth':total_depth/float(covered_bases), 'average_identity':identity/float(covered_bases),
							 'average_maf':maf/float(covered_bases)}
	# write stats
	fields = ['genome_length', 'fraction_covered', 'average_depth', 'average_identity', 'average_maf']
	outfile = open('/'.join([args['out'], 'snps_summary_stats.txt']), 'w')
	outfile.write('\t'.join(['cluster_id'] + fields)+'\n')
	for cluster_id in stats:
		record = [cluster_id] + [str(stats[cluster_id][field]) for field in fields]
		outfile.write('\t'.join(record)+'\n')

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
		infile = gzip.open(inpath)
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

def compute_phyeco_cov(args, genome_clusters, ref_to_cov, ref_to_cluster):
	""" Count number of bp mapped to each PhyEco marker gene """
	# read in set of phyeco markers for normalization
	phyeco_ids = set([])
	for line in open(args['pid_cutoffs']):
		phyeco_id, pid = line.rstrip().split()
		phyeco_ids.add(phyeco_id)
	# read in map of gene to phyeco marker
	ref_to_phyeco = {}
	for cluster_id in genome_clusters:
		inpath = '/'.join([args['db_dir'], cluster_id, 'universal_genes.txt.gz'])
		infile = gzip.open(inpath)
		next(infile)
		for line in infile:
			gene_id, phyeco_id = line.rstrip().split()
			ref_to_phyeco[gene_id] = phyeco_id
	# init phyeco coverage
	cluster_to_phyeco_to_cov = {}
	for cluster_id in genome_clusters:
		cluster_to_phyeco_to_cov[cluster_id] = {}
		for phyeco_id in phyeco_ids:
			cluster_to_phyeco_to_cov[cluster_id][phyeco_id] = 0.0
	# compute phyeco coverages
	for ref_id, phyeco_id in ref_to_phyeco.items():
		cluster_id = ref_to_cluster[ref_id]
		if phyeco_id in phyeco_ids and ref_id in ref_to_cov:
			cluster_to_phyeco_to_cov[cluster_id][phyeco_id] += ref_to_cov[ref_id]
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
		outfiles[cluster_id] = gzip.open('/'.join([outdir, '%s.cov.gz' % cluster_id]), 'w')
		outfiles[cluster_id].write('\t'.join(['ref_id', 'raw_coverage', 'normalized_coverage'])+'\n')
	# parse bam into cov files for each cluster_id
	ref_to_cov = count_mapped_bp(args)
	# compute normalization factor
	cluster_to_norm = compute_phyeco_cov(args, genome_clusters, ref_to_cov, ref_to_cluster)
	# aggregate coverages
	if args['pangenome_pid'] != '99':
		ref_to_cov = aggregate_coverage(args, ref_to_cov, genome_clusters)
	# write to output files
	for ref_id, cov in ref_to_cov.items():
		cluster_id = ref_to_cluster[ref_id]
		outfile = outfiles[cluster_id]
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
		infile = gzip.open(inpath)
		for r in Bio.SeqIO.parse(infile, 'fasta'):
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
	infile = gzip.open(inpath)
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
		infile = gzip.open(inpath)
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
	if args['snps_pileup']:
		start = time()
		if args['verbose']: print("\nRunning mpileup")
		pileup(args)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# 3d. Use mpileup to identify SNPs
	if args['snps_call']:
		start = time()
		if args['verbose']: print("\nFormatting output")
		split_vcf(args)
		format_vcf(args)
		snps_summary(args)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

