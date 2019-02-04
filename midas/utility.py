#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import io, os, stat, sys, resource, gzip, platform, bz2, Bio.SeqIO

__version__ = '1.3.0'

def which(program):
	""" Mimics unix 'which' function """
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			path = path.strip('"')
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file
	return None

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
		assert: len(batches) == threads
	"""
	import resource
	import math
	max_open = int(0.8 * resource.getrlimit(resource.RLIMIT_NOFILE)[0]) # max open files on system
	max_size = math.floor(max_open/threads) # max batch size to avoid exceeding max_open
	min_size = math.ceil(len(samples)/float(threads)) # min batch size to use all threads
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

def parallel_old(function, list, threads):
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

def parallel(function, argument_list, threads, no_results=False):
	""" Based on: https://gist.github.com/admackin/003dd646e5fadee8b8d6 """
	import multiprocessing as mp
	import signal
	import time

	def init_worker():
		signal.signal(signal.SIGINT, signal.SIG_IGN)

	pool = mp.Pool(int(threads), init_worker)

	try:
		results = []
		for arguments in argument_list:
			p = pool.apply_async(function, args=arguments)
			results.append(p)
		pool.close()

		while True:
			if all(r.ready() for r in results):
				pool.join()
				if no_results:
					return
				return [r.get() for r in results]
			time.sleep(1)

	except KeyboardInterrupt:
		pool.terminate()
		pool.join()
		sys.exit("\nKeyboardInterrupt")

def add_executables(args):
	""" Identify relative file and directory paths """
	src_dir = os.path.dirname(os.path.abspath(__file__))
	main_dir = os.path.dirname(src_dir)
	args['stream_seqs'] = '/'.join([src_dir, 'run', 'stream_seqs.py'])
	args['hs-blastn'] = '/'.join([main_dir, 'bin', platform.system(), 'hs-blastn'])
	args['bowtie2-build'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2-build'])
	args['bowtie2'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2'])
	args['samtools'] = '/'.join([main_dir, 'bin', platform.system(), 'samtools'])

	for arg in ['hs-blastn', 'stream_seqs', 'bowtie2-build', 'bowtie2', 'samtools']:
		if not os.path.isfile(args[arg]):
			sys.exit("\nError: File not found: %s\n" % args[arg])

	for arg in ['hs-blastn', 'bowtie2-build', 'bowtie2', 'samtools']:
		if not os.access(args[arg], os.X_OK):
			sys.exit("\nError: File not executable: %s\n" % args[arg])

	import subprocess as sp

	process = sp.Popen("%s view" % args['samtools'], shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	process.wait()
	if process.returncode != 0:
		err = "\nError: could not execute samtools binary: %s\n" % args['samtools']
		err += "(exited with error code %s)\n" % process.returncode
		err += "To solve this issue, follow these steps:\n"
		err += "  1) Download samtools v1.4: https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2\n"
		err += "  2) Unpack and compile the software on your system\n"
		err += "  3) Copy the new samtools binary to: %s\n" % os.path.dirname(args['samtools'])
		sys.exit(err)

	process = sp.Popen("%s -h" % args['bowtie2'], shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	process.wait()
	if process.returncode != 0:
		err = "\nError: could not execute bowtie2 binary: %s\n" % args['bowtie2']
		err += "(exited with error code %s)\n" % process.returncode
		err += "To solve this issue, follow these steps:\n"
		err += "  1) Go to https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2\n"
		err += "  2) Download bowtie2-2.3.2-linux-x86_64.zip\n"
		err += "  3) Unpack the software on your system\n"
		err += "  4) Copy the new bowtie2 binaries to: %s\n" % os.path.dirname(args['bowtie2'])
		sys.exit(err)

def auto_detect_file_type(inpath):
	""" Detect file type [fasta or fastq] of <p_reads> """
	infile = iopen(inpath)
	for line in infile:
		if line[0] == '>': return 'fasta'
		elif line[0] == '@': return 'fastq'
		else: sys.exit("Error: Filetype [fasta, fastq] of %s could not be recognized\n" % inpath)
	infile.close()

def check_compression(inpath):
	""" Check that file extension matches expected compression """
	ext = inpath.split('.')[-1]
	file = iopen(inpath)
	try:
		next(file)
		file.close()
	except:
		sys.exit("\nError: File extension '%s' does not match expected compression\n" % ext)

def check_database(args):
	if args['db'] is None:
		error = "\nError: No reference database specified\n"
		error += "Use the flag -d to specify a database,\n"
		error += "Or set the MIDAS_DB environmental variable: export MIDAS_DB=/path/to/midas/db\n"
		sys.exit(error)
	if not os.path.isdir(args['db']):
		error = "\nError: Specified reference database does not exist: %s\n" % args['db']
		error += "\nCheck that you've entered the path correctly and the database exists"
		error += "\nTo download the default database, run: MIDAS/scripts/download_ref_db.py"
		error += "\nTo build a custom database, run: MIDAS/scripts/build_midas_db.py\n"
		sys.exit(error)
	for dir in ['pangenomes', 'repgenomes']:
		path = '%s/%s' % (args['db'], dir)
		if not os.path.isdir(path):
			error = "\nError: Could not locate required database directory: %s\n" % path
			sys.exit(error)
	for file in ['metadata/species_info.tsv']:
		path = '%s/%s' % (args['db'], file)
		if not os.path.exists(path):
			error = "\nError: Could not locate required database file: %s\n" % path
			sys.exit(error)

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
	fields = next(infile).rstrip('\n').split('\t')
	for line in infile:
		values = line.rstrip('\n').split('\t')
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
		err_message = "\nError encountered executing:\n%s\n\nError message:\n%s\n" % (command, err)
		sys.exit(err_message)

def check_bamfile(args, bampath):
	""" Use samtools to check bamfile integrity """
	import subprocess as sp
	command = '%s view %s > /dev/null' % (args['samtools'], bampath)
	process = sp.Popen(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	out, err = process.communicate()
	if err.decode('ascii') != '': # need to use decode to convert to characters for python3
		err_message = "\nWarning, bamfile may be corrupt: %s\nSamtools reported this error: %s\n" % (bampath, err.rstrip())
		sys.exit(err_message)

def read_genes(species_id, db):
	""" Read in gene coordinates from features file """
	genome = read_genome(db, species_id)
	genes = []
	basename = '%s/rep_genomes/%s/genome.features' % (db, species_id)
	if os.path.exists(basename):
		fpath = basename
	elif os.path.exists(basename+'.gz'):
		fpath = basename+'.gz'
	else:
		sys.exit("\nError: rep genome for %s not found\n" % species_id)
	for gene in parse_file(fpath):
		if 'gene_type' in gene and gene['gene_type'] != 'CDS':
			continue
		else:
			gene['start'] = int(gene['start'])
			gene['end'] = int(gene['end'])
			gene['seq'] = get_gene_seq(gene, genome[gene['scaffold_id']])
			genes.append(gene)

	# sort genes
	coords = [[gene['scaffold_id'], gene['start'], -gene['end']] for gene in genes]
	indexes = sorted(range(len(coords)), key=lambda k: coords[k])
	sorted_genes = [genes[i] for i in indexes]

	return {'list':sorted_genes, 'index':0}


def read_genome(db, species_id):
	""" Read in representative genome from reference database """
	basename = '%s/rep_genomes/%s/genome.fna' % (db, species_id)
	if os.path.exists(basename):
		fpath = basename
	elif os.path.exists(basename+'.gz'):
		fpath = basename+'.gz'
	else:
		sys.exit("\nError: rep genome for %s not found\n" % species_id)
	infile = iopen(fpath)
	genome = {}
	for r in Bio.SeqIO.parse(infile, 'fasta'):
		genome[r.id] = r.seq.upper()
	infile.close()
	return genome

def get_gene_seq(gene, genome):
	""" Fetch nucleotide sequence of gene from genome """
	seq = genome[gene['start']-1:gene['end']] # 2x check this works for + and - genes
	if gene['strand'] == '-':
		return(rev_comp(seq))
	else:
		return(seq)

def complement(base):
	""" Complement nucleotide """
	d = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
	if base in d: return d[base]
	else: return base

def rev_comp(seq):
	""" Reverse complement sequence """
	return(''.join([complement(base) for base in list(seq[::-1])]))

def translate(codon):
	""" Translate individual codon """
	codontable = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
	'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
	}
	return codontable[str(codon)]

def index_replace(codon, allele, pos, strand):
	""" Replace character at index i in string x with y"""
	bases = list(codon)
	bases[pos] = allele if strand == '+' else complement(allele)
	return(''.join(bases))
