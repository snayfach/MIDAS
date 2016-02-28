#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.2'

## LIBRARIES
import argparse, sys, os, gzip, subprocess, tempfile, resource, math, shutil

## CLASSES
class Sample:
	""" Base class for samples """
	def __init__(self, dir, species_id):
		self.id = os.path.basename(dir)
		self.dir = dir
		self.average_depth = None
		self.fraction_covered = None
		self.summary_stats(species_id)

	def summary_stats(self, species_id):
		inpath = '/'.join([self.dir, 'snps_summary_stats.txt'])
		if not os.path.isfile(inpath): return
		stats = parse_summary(inpath)
		if species_id not in stats: return
		self.average_depth = stats[species_id]['average_depth']
		self.fraction_covered = stats[species_id]['fraction_covered']

class Snp:
	""" Base class for SNPs """
	def __init__(self, snpfiles, samples):
		self.data = self.store_data(snpfiles, samples)
		
	def store_data(self, snpfiles, samples):
		x = []
		for sample in samples:
			try: x.append(next(snpfiles[sample.id]))
			except StopIteration: return None
		return x

	def id(self):
		return [self.data[0]['ref_id'], self.data[0]['ref_pos']]

	def fetch(self, field, type=None):
		if type: return [type(_[field]) for _ in self.data]
		else: return [_[field] for _ in self.data]

## FUNTIONS
def parse_arguments():
	""" Parse command line arguments """
	
	parser = argparse.ArgumentParser(
		usage='%s [options]' % os.path.basename(__file__),
		description="""Merge single-nucleotide variants for an individual species across samples. Outputs include: a list of high-quality sites, an allele frequency matrix, consensus sequences for each sample, and a phylogenetic tree"""
		)
	parser.add_argument('--verbose', action='store_true', default=False,
		help="Verbosity")
	parser.add_argument('--threads', type=int, default=1,
		help="Number of threads to use")
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""input to results from 'run_phylo_cnv.py snvs'.
			see <intype> for details""")
	io.add_argument('-t', choices=['dir', 'file', 'list'], dest='intype', required=True,
		help="""input type.
			'dir': directory containing phylo_cnv results. each subdirectory should correspond to a different sample_id. for example: <directory>/<sample_id>
			'file': file containing paths to phylo_cnv results.	each line in the file should contain the full path to the results for a sample_id.
			'list': comma-separated list of paths to phylo_cnv results.
			""")
	io.add_argument('-s', dest='species_id', type=str, required=True,
		help="""species identifier. 
			a list of prevalent species can be obtained by running 'scripts/merge_species.py'.
			A map of species ids to species names can be found in 'ref_db/annotations.txt'""")
	io.add_argument('-o', dest='outdir', type=str, required=True,
		help="""output directory.
			output files: <outdir>/<species_id>.sample_ids, <outdir>/<species_id>.hq_snps, <outdir>/<species_id>.ref_freq, <outdir>/<species_id>.depth, <outdir>/<species_id>.fasta, <outdir>/<species_id>.tree""")
	#io.add_argument('-m', '--matrix', type=str,  help='reference SNP matrix')
	
	pipe = parser.add_argument_group('Pipeline options (choose one or more; default=all)')
	pipe.add_argument('--snps', default=False, action='store_true', help='identify and store list of hq snps')
	pipe.add_argument('--freq', default=False, action='store_true', help='build allele frequency & depth matrixes')
	pipe.add_argument('--cons', default=False, action='store_true', help='generate fasta file of consensus sequences')
	pipe.add_argument('--tree', default=False, action='store_true', help='build phylogenetic tree')
	
	sample = parser.add_argument_group('Sample filters\n(determine which samples are included in output)')
	sample.add_argument('--sample_depth', dest='sample_depth', type=float, default=5.0,
		help="""minimum average read depth per sample (5.0)""")
	sample.add_argument('--fract_cov', dest='fract_cov', type=float, default=0.4,
		help="""fraction of reference sites covered by at least 1 read (0.4)""")
	sample.add_argument('--max_samples', type=int,
		help="""maximum number of samples to process. useful for quick tests (use all)""")
				
	snps = parser.add_argument_group('Site filters\n(determine which reference-genome positions are included in output)')
	snps.add_argument('--site_depth', dest='site_depth', type=int, default=3,
		help="""minimum number of mapped reads per site. a high value like 20 will result in accurate allele frequencies, but may discard many sites. a low value like 1 will retain many sites but may not result in accurate allele frequencies (3)""")
	snps.add_argument('--site_prev', dest='site_prev', type=float, default=0.95,
		help="""site has at least <site_depth> coverage in at least <site_prev> proportion of samples.
			a value of 1.0 will select sites that have sufficent coverage in all samples.
			a value of 0.0 will select all sites, including those with low coverage in many samples  (0.95)""")
	snps.add_argument('--max_sites', dest='max_sites', type=int, default=float('Inf'),
		help="""maximum number of sites to include in output. useful for quick tests (use all)""")
					
	args = vars(parser.parse_args())
	args['db'] = '%s/ref_db/genome_clusters' % os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	check_args(args)
	return args

def check_args(args):
	if not os.path.isdir(args['outdir']):
		os.mkdir(args['outdir'])
	if not any([args['snps'], args['freq'], args['cons'], args['tree']]):
		args['snps'] = True
		args['freq'] = True
		args['cons'] = True
		args['tree'] = True

def print_copyright():
	print ("")
	print ("PhyloCNV: species abundance and strain-level genomic variation from metagenomes")
	print ("version %s; github.com/snayfach/PhyloCNV" % __version__)
	print ("Copyright (C) 2015 Stephen Nayfach")
	print ("Freely distributed under the GNU General Public License (GPLv3)")
	print ("")

def print_arguments(args):
	print ("-------------------------------------------------------")
	print ("Merge SNVs Parameters")
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Species identifier: %s" % args['species_id'])
	print ("Pipeline options:")
	if args['snps']: print ("  identify sites")
	if args['freq']: print ("  build allele frequency matrix")
	if args['cons']: print ("  call consensus sequences")
	if args['tree']: print ("  build phylogenetic tree")
	print ("Sample selection criteria:")
	if args['sample_depth']:
		print ("  keep samples with >=%s average coverage across reference genome" % args['sample_depth'])
	if args['fract_cov']:
		print ("  keep samples where >=%s percent of reference genome has non-zero coverage" % (100*args['fract_cov']))
	if args['max_samples']:
		print ("  analyze up to %s samples" % args['max_samples'])
	print ("Site selection criteria:")
	print ("  site must be covered by at least %s reads across %s percent of samples" % (args['site_depth'], 100*args['site_prev']))
	if args['max_sites'] != float('Inf'):
		print ("  analyze up to %s sites" % (args['max_sites']))
	print ("-------------------------------------------------------")
	print ("")

def list_samples(input, intype):
	""" Get full path to sample directories """
	if intype == 'dir':
		if not os.path.isdir(input):
			sys.exit("\nSpecified input directory does not exist:\n%s" % input)
		else:
			return([os.path.join(input, _) for _ in os.listdir(input)])
	elif intype == 'file':
		if not os.path.isfile(input):
			sys.exit("\nSpecified input file does not exist:\n%s" % input)
		else:
			return([x.rstrip() for x in open(input).readlines()])
	elif intype == 'list':
		return(input.split(','))

def parse_summary(inpath):
	""" Read in summary snp statistics for genome-clusters """
	snps_summary = {}
	infile = open(inpath)
	fields = next(infile).rstrip().split()
	for line in open(inpath):
		values = line.rstrip().split()
		rec = dict([(i,j) for i,j in zip(fields, values)])
		snps_summary[rec['cluster_id']] = rec
	return snps_summary

def write_samples(args, samples):
	""" Write sample list to file """
	outfile = open('%s/%s.sample_ids' % (args['outdir'], args['species_id']), 'w')
	outfile.write('\t'.join(['sample_id', 'average_depth', 'fraction_covered'])+'\n')
	for sample in samples:
		outfile.write(sample.id+'\n')
		outfile.write(str(sample.average_depth)+'\n')
		outfile.write(str(sample.fraction_covered)+'\n')

def identify_samples(args):
	""" Identify samples with sufficient depth for target genome-cluster """
	samples = []
	for dir in list_samples(args['input'], args['intype']):
		sample = Sample(dir, args['species_id'])
		if not sample.average_depth:
			continue
		elif float(sample.average_depth) < args['sample_depth']:
			continue
		elif float(sample.fraction_covered) < args['fract_cov']:
			continue
		else:
			samples.append(sample)
		if args['max_samples'] and len(samples) >= args['max_samples']:
			break
	if len(samples) == 0:
		sys.exit("Error: species_id failed to pass quality-control in all samples")
	else:
		write_samples(args, samples)
		print("  %s samples with species" % len(samples))
		return samples

def parse_snps(inpath):
	""" Yields formatted records from SNPs output """
	infile = gzip.open(inpath)
	next(infile)
	fields = ['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
			  'count_alleles', 'count_ref', 'count_alt', 'depth', 'ref_freq']
	for line in infile:
		yield dict([(i,j) for i,j in zip(fields, line.rstrip().split())])

def open_infiles(species_id, samples):
	""" Open SNP files for species across samples """
	infiles = {}
	for sample in samples:
		inpath = '%s/snps/%s.snps.gz' % (sample.dir, species_id)
		infiles[sample.id] = parse_snps(inpath)
	return infiles

def read_hq_snps(outdir, species_id):
	""" Read in list of HQ snps from file"""
	snps = []
	for line in open('%s/%s.hq_snps' % (outdir, species_id)):
		snps.append(line.rstrip().split()[0:2])
	return snps

def prevalence(snp, depth, freq=True):
	""" Count number/fraction of samples where SNP is HQ """
	qc = []
	for s in snp.data:
		if float(s['depth']) < depth: qc.append(False)
		else: qc.append(True)
	if freq: return float(sum(qc))/len(qc)
	else: return sum(qc)

def open_outfiles(outdir, species_id, sample_ids, index=None):
	""" Open matrices and write headers """
	outfiles = {}
	for type in ['ref_freq', 'depth', 'alt_allele']:
		if index is None: outpath = '%s/%s.%s' % (outdir, species_id, type)
		else: outpath = '%s/%s.%s.%s' % (outdir, species_id, index, type)
		outfiles[type] = open(outpath, 'w')
		outfiles[type].write('\t'.join(['ref_id','ref_pos']+sample_ids)+'\n')
	return outfiles

def batch_samples(samples, threads):
	""" Split up samples into batches
		assert: batch_size * threads < max_open
		assert: batch_size uses all threads
	"""
	import resource
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

def identify_snps(args, samples):
	""" 
	Write list of high-quality snps to disk
	  -Split up samples into batches
	  -In parallel, write temp snp lists with sample counts
	  -Loop over temp snp lists and id high quality snps
	"""
	# write temporary snp lists
	list = []
	tempdir = tempfile.mkdtemp(dir=args['outdir'])
	batches = batch_samples(samples, args['threads'])
	for index, batch in enumerate(batches):
		a = {}
		a['index'] = index
		a['samples'] = batch
		a['species_id'] = args['species_id']
		a['outdir'] = tempdir
		a['site_depth'] = args['site_depth']
		a['max_sites'] = args['max_sites']
		list.append(a)
	parallel(write_prevalence, list, args['threads'])
	# open temporary snp lists
	infiles = []
	for index, batch in enumerate(batches):
		inpath = '%s/%s.%s.hq_snps' % (tempdir, args['species_id'], index)
		infiles.append(open(inpath))
	# identify hq snps
	outfile = open('%s/%s.hq_snps' % (args['outdir'], args['species_id']), 'w')
	min_count = int(len(samples) * args['site_prev'])
	hq_count = 0
	while True:
		ref_id, ref_pos = None, None
		counts = 0
		for infile in infiles:
			try:
				ref_id, ref_pos, count = infile.next().rstrip().split()
				counts += int(count)
			except StopIteration:
				break
		if not ref_id:
			break
		elif counts >= min_count:
			outfile.write('%s\t%s\n' % (ref_id, ref_pos))
			hq_count += 1
	print "  %s high-quality snps found" % hq_count
	# clean up temporary snp lists
	shutil.rmtree(tempdir)

def write_prevalence(index, species_id, samples, site_depth, max_sites, outdir):
	""" Count number of samples that have sufficient depth across sites """
	snpfiles = open_infiles(species_id, samples)
	outpath = '%s/%s.%s.hq_snps' % (outdir, species_id, index)
	outfile = open(outpath, 'w')
	count = 0
	while count <= max_sites:
		snp = Snp(snpfiles, samples)
		if snp.data is None:
			break
		else:
			snp_prev = prevalence(snp, site_depth, freq=False)
			record = [snp.id()[0], snp.id()[1], snp_prev]
			outfile.write('\t'.join([str(_) for _ in record])+'\n')
			count += 1

def build_snp_matrix(args, samples):
	""" 
	Write snp matrices: ref_freq, depth, alt_allele
	  -Split up samples into batches
	  -In parallel, write temp files with subset of samples
	  -Loop over temp snp lists and id high quality snps
	"""
	# write temporary snp matrixes
	list = []
	tempdir = tempfile.mkdtemp(dir=args['outdir'])
	batches = batch_samples(samples, args['threads'])
	hq_snps = read_hq_snps(args['outdir'], args['species_id'])
	for index, batch in enumerate(batches):
		a = {}
		a['index'] = index
		a['samples'] = batch
		a['species_id'] = args['species_id']
		a['hq_snps'] = hq_snps
		a['outdir'] = tempdir
		list.append(a)
	parallel(build_mini_matrix, list, args['threads'])
	# open temporary snp matrixes
	infiles = {}
	for type in ['ref_freq', 'depth', 'alt_allele']:
		files = []
		for index, batch in enumerate(batches):
			inpath = '%s/%s.%s.%s' % (tempdir, args['species_id'], index, type)
			files.append(open(inpath))
		infiles[type] = files
	# merge temporary files
	sample_ids = [s.id for s in samples]
	outfiles = open_outfiles(args['outdir'], args['species_id'], sample_ids)
	for type in ['ref_freq', 'depth', 'alt_allele']:
		files = infiles[type]
		for file in files: # skip header
			next(file)
		while True:
			values = []
			try:
				for index, file in enumerate(files):
					v = next(file).rstrip().split()
					if index == 0: values += v
					else: values += v[2:]
				outfiles[type].write('\t'.join(values)+'\n')
			except StopIteration:
				break
	# remove temporary files
	shutil.rmtree(tempdir)

def build_mini_matrix(outdir, species_id, samples, hq_snps, index):
	""" Build SNP matrices using a subset of total samples """
	sample_ids = [s.id for s in samples]
	outfiles = open_outfiles(outdir, species_id, sample_ids, index)
	snpfiles = open_infiles(species_id, samples)
	nsnps = len(hq_snps)
	pos = 0 # index position in hq_snps
	while True:
		snp = Snp(snpfiles, samples) # snp across samples
		if snp.data is None: # eof
			break
		elif pos == nsnps: # no more snps
			break
		elif snp.id() != hq_snps[pos]: # snp not in list
			continue
		else: # snp in list
			pos += 1
			id = snp.id()
			for field in ['ref_freq', 'depth', 'alt_allele']:
				values = snp.fetch(field, str)
				outfiles[field].write('\t'.join(id)+'\t'+'\t'.join(values)+'\n')

def write_consensus(args, samples):
	""" Write consensus sequences from samples """
	outfile = open('%s/%s.fasta' % (args['outdir'], args['species_id']), 'w')
	hq_snps = read_hq_snps(args['outdir'], args['species_id']) # get snp list
	n = len(hq_snps)
	for sample in samples:
		i = 0 # index position in hq_snps
		outfile.write('>%s\n' % sample.id) # sequence header
		inpath = '%s/snps/%s.snps.gz' % (sample.dir, args['species_id'])
		for snp in parse_snps(inpath):
			if i == n: # no more snps in list
				break
			elif [snp['ref_id'], snp['ref_pos']] != hq_snps[i]: # snp not in list
				continue
			else: # snp in list
				outfile.write('-' if snp['cons_allele'] == 'NA' else snp['cons_allele']) # write base
				i += 1
		outfile.write('\n')

def build_tree(args):
	"""	Use FastTree to build phylogenetic tree of consensus sequences """
	inpath = '%s/%s.fasta' % (args['outdir'], args['species_id'])
	outpath = '%s/%s.tree' % (args['outdir'], args['species_id'])
	p = subprocess.Popen('FastTree -nt -boot 100 < %s > %s' % (inpath, outpath), shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	out, err = p.communicate()


## MAIN
if __name__ == '__main__':

	args = parse_arguments()
	print_copyright()
	print_arguments(args)

	print("Identifying hq samples with species")
	samples = identify_samples(args)
	
	if args['snps']:
		print("Identifying and writing hq snps")
		identify_snps(args, samples)

	if args['freq']:
		print("Writing allele frequencies & depths")
		build_snp_matrix(args, samples)

	if args['cons']:
		print("Writing consensus sequences")
		write_consensus(args, samples)

	if args['tree']:
		print("Building phylogenetic tree")
		build_tree(args)

			
