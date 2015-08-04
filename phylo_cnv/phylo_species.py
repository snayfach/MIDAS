#!/usr/bin/python

# Libraries
# ---------
import sys
import os
import gzip
import subprocess
import Bio.SeqIO
import random
import numpy as np
import time
from tempfile import mkstemp

# Functions
# ---------
def parse_relative_paths(args):
	""" Identify relative file and directory paths """
	paths = {}
	paths['main'] = os.path.dirname(os.path.abspath(__file__))
	paths['blastn'] = '/'.join([paths['main'],'bin','blastn'])
	assert(os.path.isfile(paths['blastn']))
	paths['cluster_ids'] = '/'.join([paths['main'],'data','cluster_annotations.txt'])
	assert(os.path.isfile(paths['cluster_ids']))
	paths['gene_length'] = '/'.join([paths['main'],'data','gene_length.txt'])
	assert(os.path.isfile(paths['gene_length']))
	paths['marker_cutoffs'] = '/'.join([paths['main'],'data','pid_cutoffs.txt'])
	assert(os.path.isfile(paths['marker_cutoffs']))
	paths['tempfile'] = mkstemp(dir=args['temp_dir'])[1] if args['temp_dir'] else mkstemp()[1]
	paths['db'] = '/'.join([paths['main'],'data','phyeco.blastdb'])
	return paths

def auto_detect_file_type(inpath):
	""" Detect file type [fasta or fastq] of <p_reads> """
	for line in iopen(inpath):
		if line[0] == '>': return 'fasta'
		elif line[0] == '@': return 'fastq'
		else: sys.exit("Filetype [fasta, fastq] of %s could not be recognized" % inpath)

def auto_detect_fastq_format(inpath):
	""" Use first 50,000 reads to detect quality encoding """
	max_reads = 50000
	formats = ['fastq-illumina', 'fastq-solexa', 'fastq-sanger']
	for format in formats:
		try:
			infile = iopen(inpath)
			for index, rec in enumerate(Bio.SeqIO.parse(infile, format)):
				if index == max_reads: break
			return format
		except:
			infile.close()
	sys.exit("Could not determine FASTQ quality encoding")

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

def quality_filter(rec, quality_type, min_length, min_quality, max_n):
	""" Return true if read fails QC """
	# check for Ns
	if sum([1 if b == 'N' else 0 for b in rec.seq])/float(len(rec.seq)) > max_n:
		return True
	# check length
	if len(rec.seq) < min_length:
		return True
	# check mean quality
	quality = rec.letter_annotations[quality_type]
	if quality and np.mean(quality) < min_quality:
		return True
	# read passed QC
	return False

#def fastq_to_fasta(args, paths):
#	""" Sample high quality reads from seqfile(s) """
#	outfile = open(paths['tempfile'], 'w')
#	index = 0
#	for file in args['inpaths']:
#		if args['indir']: inpath = os.path.join(args['indir'], file)
#		else: inpath = file
#		file_format = auto_detect_file_type(inpath)
#		if file_format == 'fastq':
#			fastq_format = auto_detect_fastq_format(inpath)
#			quality_type = 'solexa_quality' if fastq_format == 'fastq-solexa' else 'phred_quality'
#		for rec in Bio.SeqIO.parse(iopen(inpath), fastq_format if file_format == 'fastq' else 'fasta'):
#			if file_format == 'fastq' and quality_filter(rec, quality_type, args['min_length'], args['min_quality'], args['max_n']):
#				continue
#			else:
#				outfile.write('>'+str(rec.id)+'\n'+str(rec.seq)+'\n')
#				index += 1
#				if index == args['nreads']: return

def fastq_to_fasta(args, paths):
	""" Sample high quality reads from seqfile(s) """
	outfile = open(paths['tempfile'], 'w')
	index = 0
	file_format = auto_detect_file_type(args['inpath'])
	if file_format == 'fastq':
		fastq_format = auto_detect_fastq_format(args['inpath'])
		quality_type = 'solexa_quality' if fastq_format == 'fastq-solexa' else 'phred_quality'
	for rec in Bio.SeqIO.parse(iopen(args['inpath']), fastq_format if file_format == 'fastq' else 'fasta'):
		if file_format == 'fastq' and quality_filter(rec, quality_type, args['min_length'], args['min_quality'], args['max_n']):
			continue
		else:
			outfile.write('>'+str(rec.id)+'\n'+str(rec.seq)+'\n')
			index += 1
			if index == args['nreads']: return

def map_reads_blast(args, paths):
	""" Use blastn to map reads in fasta file to marker database """
	command = " %(blastn)s -query %(query)s -db %(db)s -out %(out)s -outfmt 6 -num_threads %(threads)s"
	arguments = {'blastn':paths['blastn'], 'query':paths['tempfile'], 'db':paths['db'], 'out':'/dev/stdout', 'threads':args['threads']}
	process = subprocess.Popen(command % arguments, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = process.communicate()
	return stdout

def parse_blast(blastout):
	""" Yield formatted record from BLAST m8 file """
	formats = [str,str,float,int,float,float,float,float,float,float,float,float]
	fields = ['query','target','pid','aln','mis','gaps','qstart','qend','tstart','tend','evalue','score']
	for line in blastout.split('\n')[0:-1]:
		values = line.rstrip().split()
		yield dict([(field, format(value)) for field, format, value in zip(fields, formats, values)])

def find_best_hits(blastout, paths):
	""" Find top scoring alignment for each read """
	best_hits = {}
	marker_cutoffs = get_markers(paths)
	for aln in blastout:
		marker_id = aln['target'].split('_')[-1]
		if aln['pid'] < marker_cutoffs[marker_id]: # does not meet marker cutoff
			continue
		elif aln['query'] not in best_hits: # record aln
			best_hits[aln['query']] = [aln]
		elif best_hits[aln['query']][0]['score'] == aln['score']: # add aln
			best_hits[aln['query']] += [aln]
		elif best_hits[aln['query']][0]['score'] < aln['score']: # update aln
			best_hits[aln['query']] = [aln]
	return best_hits.values()

def assign_unique(args, paths, alns):
	""" Count the number of uniquely mapped reads to each genome cluster """
	unique_alns = dict([(x.rstrip().split()[0],[]) for x in open(paths['cluster_ids']).readlines()])
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
		print("\tuniquely mapped reads: %s" % unique)
		print("\tambiguously mapped reads: %s" % non_unique)
	return unique_alns

def assign_non_unique(args, paths, alns, unique_alns):
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

def get_markers(paths):
	""" Read in optimal mapping parameters for marker genes """
	marker_cutoffs = {}
	infile = open(paths['marker_cutoffs'])
	for line in infile:
		marker_id, min_pid = line.rstrip().split()
		marker_cutoffs[marker_id] = float(min_pid)
	return marker_cutoffs

def impute_missing_args(args):
	""" Fill in defaults for missing arguments """
	if 'temp_dir' not in args:
		args['temp_dir'] = None
	if 'indir' not in args:
		args['indir'] = None
	if 'max_n' not in args:
		args['max_n'] = 1.0
	if 'min_quality' not in args:
		args['min_quality'] = 0
	if 'min_length' not in args:
		args['min_length'] = 0
	if 'verbose' not in args:
		args['verbose'] = False
	if 'threads' not in args:
		args['threads'] = 1
	if 'normalize' not in args:
		args['normalize'] = False
	if 'nreads' not in args:
		args['nreads'] = float('Inf')
	return args

def estimate_mix_props(alns, paths):
	""" Count the number of uniquely mapped reads to each genome cluster """
	unique_counts = dict([(x.rstrip().split()[0],0) for x in open(paths['cluster_ids']).readlines()])
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

def read_gene_lengths(paths):
	""" Read in total gene length per cluster_id """
	total_gene_length = dict([(x.rstrip().split()[0],0) for x in open(paths['cluster_ids']).readlines()])
	for line in open(paths['gene_length']):
		cluster_id = line.split()[0].split('_')[0]
		gene_length = int(line.rstrip().split()[1])
		total_gene_length[cluster_id] += gene_length
	return total_gene_length

def normalize_counts(cluster_alns, total_gene_length, total_genomes):
	""" Normalize counts by gene length and sum contrain """
	# norm by gene length, compute rpkg, compute cov
	cluster_abundance = {}
	for cluster_id, alns in cluster_alns.items():
		cluster_abundance[cluster_id] = {}
		# compute various abundance metrics
		if len(alns) > 0:
			reads = len(alns)
			bp = sum([aln['aln'] for aln in alns])
			rpkg = float(reads)/(total_gene_length[cluster_id]/1000.0)/total_genomes
			cov = float(bp)/total_gene_length[cluster_id]
			prop_cov = cov/total_genomes
		else:
			reads, bp, rpkg, cov, prop_cov = 0.0, 0.0, 0.0, 0.0, 0.0
		# store results
		cluster_abundance[cluster_id] = {'reads':reads,'bp':bp,'rpkg':rpkg,'cov':cov,'prop_cov':prop_cov}
	# compute relative abundance
	total_rpkg = sum([_['rpkg'] for _ in cluster_abundance.values()])
	for cluster_id in cluster_abundance.keys():
		rpkg = cluster_abundance[cluster_id]['rpkg']
		cluster_abundance[cluster_id]['rel_abun'] = rpkg/total_rpkg
	return cluster_abundance

def write_abundance(outpath, cluster_abundance):
	""" Write cluster results to specified output file """
	outfile = open(outpath, 'w')
	fields = ['cluster_id', 'cluster_name', 'reads', 'bp', 'rpkg', 'cov', 'prop_cov', 'rel_abun']
	outfile.write('\t'.join(fields)+'\n')
	for cluster_id, values in cluster_abundance.items():
		reads = values['reads']
		bp = values['bp']
		rpkg = values['rpkg']
		cov = values['cov']
		prop_cov = values['prop_cov']
		rel_abun = values['rel_abun']
		record = [cluster_id, reads, bp, rpkg, cov, prop_cov, rel_abun]
		outfile.write('\t'.join([str(x) for x in record])+'\n')

def write_summary(outpath, cluster_summary):
	""" Write cluster summary to specified output file """
	outfile = open(outpath, 'w')
	fields = ['total_reads', 'total_bp', 'total_ags', 'total_coverage', 'classified_reads', 'classified_bp', 'classified_coverage', 'fraction_coverage']
	outfile.write('\t'.join(fields)+'\n')
	total_reads = cluster_summary['total_reads']
	total_bp = cluster_summary['total_bp']
	total_ags = cluster_summary['total_ags']
	total_coverage = cluster_summary['total_coverage']
	classified_reads = cluster_summary['classified_reads']
	classified_bp = cluster_summary['classified_bp']
	classified_coverage = cluster_summary['classified_coverage']
	fraction_coverage = cluster_summary['fraction_coverage']
	record = [total_reads, total_bp, total_ags, total_coverage, classified_reads, classified_bp, classified_coverage, fraction_coverage]
	outfile.write('\t'.join([str(x) for x in record])+'\n')

def count_reads_bases(inpath):
	""" Count the number of base pairs in infile """
	count_reads = 0
	count_bp = 0
	for rec in Bio.SeqIO.parse(inpath, 'fasta'):
		count_reads += 1
		count_bp += len(rec.seq)
	return (count_reads, count_bp)

def estimate_species_abundance(args):
	""" Run entire pipeline """

	# impute missing args & get relative file paths
	args = impute_missing_args(args)
	paths = parse_relative_paths(args)

	# sample hq reads
	start = time.time()
	if args['verbose']: print("Sampling HQ reads")
	fastq_to_fasta(args, paths)
	if args['verbose']: print("\t %ss" % round(time.time() - start))
	
	# estimate AGS
	if args['normalize']:
		from microbe_census import microbe_census
		start = time.time()
		if args['verbose']: print("Estimating AGS")
		ags = microbe_census.run_pipeline({'seqfile':paths['tempfile']})[0]
		if args['verbose']: print("\t %ss" % round(time.time() - start))
	else:
		ags = 3000000.0
	total_reads, total_bp = count_reads_bases(paths['tempfile'])
	total_genomes = total_bp/ags

	# align reads
	start = time.time()
	if args['verbose']: print("Aligning reads")
	blastout = parse_blast(map_reads_blast(args, paths))
	if args['verbose']: print("\t %ss" % round(time.time() - start))
	
	# find best hit for each read
	start = time.time()
	if args['verbose']: print("Classifying reads")
	best_hits = find_best_hits(blastout, paths)
	unique_alns = assign_unique(args, paths, best_hits)
	cluster_alns = assign_non_unique(args, paths, best_hits, unique_alns)
	if args['verbose']: print("\t %ss" % round(time.time() - start))
	
	# estimate genome cluster abundance
	start = time.time()
	if args['verbose']: print("Estimating cluster abundance")
	total_gene_length = read_gene_lengths(paths)
	cluster_abundance = normalize_counts(cluster_alns, total_gene_length, total_genomes)
	if args['verbose']: print("\t %ss" % round(time.time() - start))
	
	# clean up and return results
	os.remove(paths['tempfile'])
	
	# cluster summary
	cluster_summary = {'total_reads':total_reads, 'total_bp':total_bp, 'total_ags':ags, 'total_coverage':total_genomes}
	cluster_summary['classified_reads'] = sum([x['reads'] for x in cluster_abundance.values()])
	cluster_summary['classified_bp'] = sum([x['bp'] for x in cluster_abundance.values()])
	cluster_summary['classified_coverage'] = sum([x['cov'] for x in cluster_abundance.values()])
	cluster_summary['fraction_coverage'] = sum([x['prop_cov'] for x in cluster_abundance.values()])
	
	# return results
	return cluster_abundance, cluster_summary



