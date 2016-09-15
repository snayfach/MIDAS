#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import io, os, stat, sys, resource, gzip, platform, subprocess, bz2

__version__ = '1.0.0'

def print_copyright(log=None):
	lines = []
	lines.append("")
	lines.append("MIDAS: Metagenomic Intra-species Diversity Analysis System")
	lines.append("version %s; github.com/snayfach/MIDAS" % __version__)
	lines.append("Copyright (C) 2015-2016 Stephen Nayfach")
	lines.append("Freely distributed under the GNU General Public License (GPLv3)")
	lines.append("")
	if log is not None: log.write('\n'.join(lines)+'\n')
	sys.stdout.write('\n'.join(lines)+'\n')

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

def add_executables(args):
	""" Identify relative file and directory paths """
	src_dir = os.path.dirname(os.path.abspath(__file__))
	main_dir = os.path.dirname(src_dir)
	args['stream_bam'] = '/'.join([src_dir, 'run', 'stream_bam.py'])
	args['stream_seqs'] = '/'.join([src_dir, 'run', 'stream_seqs.py'])
	args['hs-blastn'] = '/'.join([main_dir, 'bin', platform.system(), 'hs-blastn'])
	args['bowtie2-build'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2-build'])
	args['bowtie2'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2'])
	args['samtools'] = '/'.join([main_dir, 'bin', platform.system(), 'samtools'])
	for arg in ['hs-blastn', 'stream_seqs', 'bowtie2-build', 'bowtie2', 'samtools', 'stream_bam']:
		if not os.path.isfile(args[arg]):
			sys.exit("File not found: %s" % args[arg])
	for arg in ['hs-blastn', 'bowtie2-build', 'bowtie2', 'samtools']:
		if not is_executable(args[arg]):
			sys.exit("File not executable: %s" % args[arg])

def read_ref_to_cluster(inpath):
	""" Read in map of scaffold id to genome-cluster id """
	ref_to_cluster = {}
	infile = open(inpath)
	for line in infile:
		ref_id, cluster_id = line.rstrip().split()
		ref_to_cluster[ref_id] = cluster_id
	infile.close()
	return ref_to_cluster

def is_executable(f):
	""" Check if file is executable by all """
	st = os.stat(f)
	return bool(st.st_mode & stat.S_IXOTH)

def auto_detect_file_type(inpath):
	""" Detect file type [fasta or fastq] of <p_reads> """
	infile = iopen(inpath)
	for line in infile:
		if line[0] == '>': return 'fasta'
		elif line[0] == '@': return 'fastq'
		else: sys.exit("Filetype [fasta, fastq] of %s could not be recognized" % inpath)
	infile.close()

def check_compression(inpath):
	""" Check that file extension matches expected compression """
	ext = inpath.split('.')[-1]
	file = iopen(inpath)
	try:
		next(file)
		file.close()
	except:
		sys.exit("\nError: File extension '%s' does not match expected compression" % ext)

def iopen(inpath, mode='r'):
	""" Open input file for reading regardless of compression [gzip, bzip] or python version """
	ext = inpath.split('.')[-1]
	# Python2
	if sys.version_info[0] == 2:
		if ext == 'gz': return gzip.open(inpath, mode)
		elif ext == 'bz2': return bz2.BZ2File(inpath, mode)
		else: return open(inpath, mode)
	# Python3
	elif sys.version_info[0] == 3:
		if ext == 'gz': return io.TextIOWrapper(gzip.open(inpath, mode))
		elif ext == 'bz2': return bz2.BZ2File(inpath, mode)
		else: return open(inpath, mode)

def parse_file(inpath):
	""" Yields records from tab-delimited file with header """
	infile = iopen(inpath)
	fields = next(infile).rstrip().split('\t')
	for line in infile:
		values = line.rstrip().split('\t')
		if len(fields) == len(values):
			yield dict([(i,j) for i,j in zip(fields, values)])
	infile.close()

def max_mem_usage():
	""" Return max mem usage (Gb) of self and child processes """
	max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
	if platform.system() == 'Linux':
		return round((max_mem_self + max_mem_child)/float(1e6), 2)
	else:
		return round((max_mem_self + max_mem_child)/float(1e9), 2)

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
	if err.decode('ascii') != '': # need to use decode to convert to characters for python3
		err_message = "\nWarning, bamfile may be corrupt: %s\nSamtools reported this error: %s\n" % (bampath, err.rstrip())
		sys.exit(err_message)
