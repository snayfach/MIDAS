import os
import stat
import sys
import resource
import gzip
import platform
import subprocess

def is_executable(f):
	""" Check if file is executable by all """
	st = os.stat(f)
	return bool(st.st_mode & stat.S_IXOTH)

def add_executables(args):
	""" Identify relative file and directory paths """
	main_dir = os.path.dirname(os.path.abspath(__file__))
	args['hs-blastn'] = '/'.join([main_dir, 'bin', platform.system(), 'hs-blastn'])
	args['fa_to_fq'] = '/'.join([main_dir, 'fa_to_fq.py'])
	args['bowtie2-build'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2-build'])
	args['bowtie2'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2'])
	args['samtools'] = '/'.join([main_dir, 'bin', platform.system(), 'samtools'])
	args['filter_bam'] = '/'.join([main_dir, 'filter_bam.py'])
	for arg in ['hs-blastn', 'fa_to_fq', 'bowtie2-build', 'bowtie2', 'samtools', 'filter_bam']:
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

def parse_file(inpath):
	""" Yields formatted records from SNPs output """
	infile = gzip.open(inpath)
	fields = next(infile).rstrip().split()
	for line in infile:
		yield dict([(i,j) for i,j in zip(fields, line.rstrip().split())])

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
	args['filter_bam'] = '/'.join([main_dir, 'filter_bam.py'])

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
