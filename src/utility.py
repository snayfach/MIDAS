#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, stat, sys, resource, gzip, platform, subprocess

__version__ = '1.0.0'

def print_copyright():
	print ("")
	print ("MIDAS: Metagenomic Intra-species Diversity Analysis System")
	print ("version %s; github.com/snayfach/PhyloCNV" % __version__)
	print ("Copyright (C) 2015-2016 Stephen Nayfach")
	print ("Freely distributed under the GNU General Public License (GPLv3)")
	print ("")

def batch_samples(samples, threads):
	""" Split up samples into batches
		assert: batch_size * threads < max_open
		assert: batch_size uses all threads
	"""
	import resource
	import math
	max_open = int(0.8 * resource.getrlimit(resource.RLIMIT_NOFILE)[0]) # max open files on system
	max_size = math.floor(max_open/threads) # max batch size to avoid exceeding max_open
	min_size = math.ceil(len(samples)/threads) # min batch size to use all threads
	size = min(min_size, max_size)
	batches = []
	batch = []
	for sample in samples:
		batch.append(sample)
		if len(batch) >= size:
			batches.append(batch)
			batch = []
	if len(batch) > 0: batches.append(batch)
	return batches

def parallel(function, list, threads):
	""" Run function using multiple threads """
	from multiprocessing import Process
	from time import sleep
	processes = []
	for pargs in list: # run function for each set of args in args_list
		p = Process(target=function, kwargs=pargs)
		processes.append(p)
		p.start()
		while len(processes) >= threads: # control number of active processes
			sleep(1)
			indexes = []
			for index, process in enumerate(processes): # keep alive processes
				if process.is_alive(): indexes.append(index)
			processes = [processes[i] for i in indexes]
	while len(processes) > 0: # wait until no active processes
		sleep(1)
		indexes = []
		for index, process in enumerate(processes):
			if process.is_alive(): indexes.append(index)
		processes = [processes[i] for i in indexes]

def add_ref_db(args):
	""" Add path to reference database """
	if 'db' not in args or not args['db']:
		script_path = os.path.abspath(__file__)
		script_dir = os.path.dirname(script_path)
		main_dir = os.path.dirname(script_dir)
		args['db'] = '%s/ref_db' % main_dir
	if not os.path.isdir(args['db']):
		sys.exit("Could not locate reference database: %s" % args['db'])
	return args

def read_ref_to_cluster(args, type):
	""" Read in map of scaffold id to genome-cluster id """
	ref_to_cluster = {}
	dir = 'snps' if type == 'genomes' else 'genes'
	infile = open('/'.join([args['outdir'], '%s/db/%s.map' % (dir, type)]))
	for line in infile:
		ref_id, cluster_id = line.rstrip().split()
		ref_to_cluster[ref_id] = cluster_id
	infile.close()
	return ref_to_cluster

def is_executable(f):
	""" Check if file is executable by all """
	st = os.stat(f)
	return bool(st.st_mode & stat.S_IXOTH)

def add_executables(args):
	""" Identify relative file and directory paths """
	main_dir = os.path.dirname(os.path.abspath(__file__))
	args['hs-blastn'] = '/'.join([main_dir, 'bin', platform.system(), 'hs-blastn'])
	args['stream_seqs'] = '/'.join([main_dir, 'stream_seqs.py'])
	args['bowtie2-build'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2-build'])
	args['bowtie2'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2'])
	args['samtools'] = '/'.join([main_dir, 'bin', platform.system(), 'samtools'])
	args['stream_bam'] = '/'.join([main_dir, 'stream_bam.py'])
	for arg in ['hs-blastn', 'stream_seqs', 'bowtie2-build', 'bowtie2', 'samtools', 'stream_bam']:
		if not os.path.isfile(args[arg]):
			sys.exit("File not found: %s" % args[arg])
	for arg in ['hs-blastn', 'bowtie2-build', 'bowtie2', 'samtools']:
		if not is_executable(args[arg]):
			sys.exit("File not executable: %s" % args[arg])

def add_data_files(args):
	""" Identify relative file and directory paths """
	main_dir = os.path.dirname(os.path.abspath(__file__))
	args['cluster_ids'] = '/'.join([main_dir,'data/cluster_annotations.txt'])
	args['gene_length'] = '/'.join([main_dir,'data/gene_length.txt'])
	args['pid_cutoffs'] = '/'.join([main_dir,'data/pid_cutoffs.txt'])
	args['bad_gcs'] = '/'.join([main_dir,'data/bad_cluster_ids.txt'])
	for arg in ['cluster_ids', 'gene_length', 'pid_cutoffs', 'bad_gcs']:
		if not os.path.isfile(args[arg]):
			sys.exit("File not found: %s" % args[arg])

def auto_detect_file_type(inpath):
	""" Detect file type [fasta or fastq] of <p_reads> """
	infile = iopen(inpath)
	for line in infile:
		if line[0] == '>': return 'fasta'
		elif line[0] == '@': return 'fastq'
		else: sys.exit("Filetype [fasta, fastq] of %s could not be recognized" % inpath)
	infile.close()

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

def parse_file(inpath):
	""" Yields records from tab-delimited file with header """
	ext = inpath.split('.')[-1]
	infile = gzip.open(inpath) if ext=='gz' else open(inpath)
	fields = next(infile).rstrip().split()
	for line in infile:
		yield dict([(i,j) for i,j in zip(fields, line.rstrip().split())])
	infile.close()

def max_mem_usage():
	""" Return max mem usage (Gb) of self and child processes """
	max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
	if platform.system() == 'Linux':
		return round((max_mem_self + max_mem_child)/float(1e6), 2)
	else:
		return round((max_mem_self + max_mem_child)/float(1e9), 2)

def add_paths(args):
	""" Add paths to external files and binaries """
	main_dir = os.path.dirname(os.path.abspath(__file__))
	args['bowtie2-build'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2-build'])
	args['bowtie2'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2'])
	args['samtools'] = '/'.join([main_dir, 'bin', platform.system(), 'samtools'])
	args['pid_cutoffs'] = '/'.join([main_dir, 'data', 'pid_cutoffs.txt'])
	args['bad_gcs'] = '/'.join([main_dir, 'data', 'bad_cluster_ids.txt'])
	args['stream_bam'] = '/'.join([main_dir, 'stream_bam.py'])

def check_exit_code(process, command):
	""" Capture stdout, stderr. Check unix exit code and exit if non-zero """
	out, err = process.communicate()
	if process.returncode != 0:
		err_message = "\nError encountered executing:\n%s\n\nError message:\n%s" % (command, err)
		sys.exit(err_message)

def check_bamfile(args, bampath):
	""" Use samtools to check bamfile integrity """
	command = '%s view %s > /dev/null' % (args['samtools'], bampath)
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	if err != '':
		err_message = "\nWarning, bamfile may be corrupt: %s\nSamtools reported this error: %s\nTry rerunning with --align" % (bampath, err.rstrip())
		sys.exit(err_message)
