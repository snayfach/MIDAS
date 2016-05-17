#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

## LIBRARIES
import argparse, sys, os, gzip, subprocess, tempfile, shutil
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import annotate_snps, merge, utility

## CLASSES
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

	def ref_allele(self):
		return self.data[0]['ref_allele']
		
	def fetch(self, field, type=None):
		if type: return [type(_[field]) for _ in self.data]
		else: return [_[field] for _ in self.data]

## FUNTIONS
def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: merge_snps.py [options]

Description: merge single-nucleotide variants for an individual species across samples
Input: list of sample directories
Output: core-genome SNPs, SNP annotations, SNP allele frequency matrix, SNP alternate alleles, SNP depth,
        core-genome consensus sequences, and a phylogenetic tree
""",
		epilog="""Examples:
1) Merge results for species 57955. Provide list of paths to sample directories:
merge_snps.py -s 57955 -o outdir/57955 -i sample_1,sample_2 -t list

2) Only use samples with >15x average depth and only use sites covered by >=10 reads in at least >=95% of samples:
merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --sample_depth 15 --site_depth 10 --site_prev 95

3) Just identify core-genome sites and build matrixes; do not build consensus seqs or phylogenetic trees:
merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --snps --freq

4) Run a quick test:
merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --max_samples 10 --max_sites 1000
""")
	
	parser.add_argument('--threads', type=int, default=1,
		help="Number of threads to use")
	
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""input to sample directories output by run_phylo_cnv.py genes
see '-t' for details""")
	io.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True,
		help="""'list': -i is a comma-separated list of paths to sample directories (ex: /sample1,/sample2)
'dir': -i is a  directory containing all samples (ex: /samples_dir)
'file': -i is a file containing paths to sample directories (ex: sample_paths.txt)
""")
	io.add_argument('-s', dest='species_id', type=str, required=True,
		help="""species identifier
a list of prevalent species can be obtained by running 'merge_species.py'
a map of species ids to species names can be found in 'ref_db/annotations.txt'""")
	io.add_argument('-o', type=str, dest='outdir', required=True,
		help="""output directory""")
			
	pipe = parser.add_argument_group("Pipeline options (choose one or more; default=all)")
	pipe.add_argument('--snps', default=False, action='store_true',
		help="identify and store list of hq snps")
	pipe.add_argument('--freq', default=False, action='store_true',
		help="build allele frequency & depth matrixes")
	pipe.add_argument('--cons', default=False, action='store_true',
		help="generate fasta file of consensus sequences")
	pipe.add_argument('--tree', default=False, action='store_true',
		help="build phylogenetic tree")
	
	sample = parser.add_argument_group("Sample filters (select subset of samples from INPUT)")
	sample.add_argument('--sample_depth', dest='sample_depth', type=float, default=5.0,
		help="""minimum average read depth per sample (5.0)""")
	sample.add_argument('--fract_cov', dest='fract_cov', type=float, default=0.4,
		help="""fraction of reference sites covered by at least 1 read (0.4)""")
	sample.add_argument('--max_samples', type=int,
		help="""maximum number of samples to process.
useful for quick tests (use all)""")
	
	snps = parser.add_argument_group("Site filters (select subset of genomic sites from INPUT)")
	snps.add_argument('--site_depth', type=int, default=3,
		help="""minimum number of mapped reads per site.
a high value like 20 will result in accurate allele frequencies, but may discard many sites.
a low value like 1 will retain many sites but may not result in accurate allele frequencies (3)""")
	snps.add_argument('--site_prev', type=float, default=0.95,
		help="""site has at least <site_depth> coverage in at least <site_prev> proportion of samples.
a value of 1.0 will select sites that have sufficent coverage in all samples.
a value of 0.0 will select all sites, including those with low coverage in many samples 
NAs recorded for included sites with less than <site_depth> in a sample (0.95)""")
	snps.add_argument('--max_sites', type=int, default=float('Inf'),
		help="""maximum number of sites to include in output.
useful for quick tests (use all)""")
	args = vars(parser.parse_args())
	args = utility.add_ref_db(args)
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

def print_arguments(args):
	print ("===========Parameters===========")
	print ("Script: merge_snps.py")
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
	print ("")

def write_samples(args, samples):
	""" Write sample list to file """
	outfile = open('%s/%s.snps.sample_ids' % (args['outdir'], args['species_id']), 'w')
	outfile.write('\t'.join(['sample_id', 'average_depth', 'fraction_covered'])+'\n')
	for sample in samples:
		outfile.write(sample.id+'\n')
		outfile.write(str(sample.average_depth)+'\n')
		outfile.write(str(sample.fraction_covered)+'\n')

def identify_samples(args):
	""" Identify samples with sufficient depth for target genome-cluster """
	samples = []
	for sample in merge.load_samples(args):
		if not sample.paths['snps']:
			sys.stderr.write("Warning: no snps output for sample: %s\n" % sample.dir)
			continue
		stats = merge.read_stats(sample.paths['snps'])
		if args['species_id'] not in stats:
			continue
		elif float(stats[args['species_id']]['average_depth']) < args['sample_depth']:
			continue
		elif float(stats[args['species_id']]['fraction_covered']) < args['fract_cov']:
			continue
		else:
			samples.append(sample)
		if args['max_samples'] and len(samples) >= args['max_samples']:
			break
	if len(samples) == 0:
		sys.exit("\nError: no samples with sufficient coverage for species_id.\nTry running with more lenient parameters")
	else:
#		write_samples(args, samples)
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
		inpath = '%s/snps/snps/%s.snps.gz' % (sample.dir, species_id)
		infiles[sample.id] = parse_snps(inpath)
	return infiles

def read_hq_snps(outdir, species_id):
	""" Read in list of HQ snps from file"""
	snps = []
	infile = open('%s/%s.snps.list' % (outdir, species_id))
	next(infile)
	for line in infile:
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
		else: outpath = '%s/%s.snps.%s.%s' % (outdir, species_id, type, index)
		outfiles[type] = open(outpath, 'w')
		outfiles[type].write('\t'.join(['ref_id','ref_pos']+sample_ids)+'\n')
	return outfiles

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
	batches = utility.batch_samples(samples, args['threads'])
	for index, batch in enumerate(batches):
		a = {}
		a['index'] = index
		a['samples'] = batch
		a['species_id'] = args['species_id']
		a['outdir'] = tempdir
		a['site_depth'] = args['site_depth']
		a['max_sites'] = args['max_sites']
		list.append(a)
	utility.parallel(write_prevalence, list, args['threads'])
	# open temporary snp lists
	infiles = []
	for index, batch in enumerate(batches):
		infile = open('%s/%s.snps.list.%s' % (tempdir, args['species_id'], index))
		next(infile)
		infiles.append(infile)
	# identify hq snps
	outfile = open('%s/%s.snps.list' % (args['outdir'], args['species_id']), 'w')
	header = ['ref_id', 'ref_pos', 'ref_allele', 'count_samples']
	outfile.write('\t'.join(header)+'\n')
	min_samples = int(len(samples) * args['site_prev'])
	count_sites = 0
	while True:
		ref_id, ref_pos, ref_allele = None, None, None
		count_samples = 0
		for infile in infiles:
			try:
				ref_id, ref_pos, ref_allele, count = infile.next().rstrip().split()
				count_samples += int(count)
			except StopIteration:
				break
		if not ref_id:
			break
		elif count_samples >= min_samples:
			count_sites += 1
			values = [ref_id, ref_pos, ref_allele, count_samples]
			outfile.write('\t'.join([str(_) for _ in values])+'\n')
	print "  %s high-quality snps found" % count_sites
	# clean up temporary snp lists
	for file in infiles: file.close()
	shutil.rmtree(tempdir, ignore_errors=True)

def write_prevalence(index, species_id, samples, site_depth, max_sites, outdir):
	""" Count number of samples that have sufficient depth at each site """
	snpfiles = open_infiles(species_id, samples)
	outpath = '%s/%s.snps.list.%s' % (outdir, species_id, index)
	outfile = open(outpath, 'w')
	header = ['ref_id', 'ref_pos', 'ref_allele', 'count_samples']
	outfile.write('\t'.join(header)+'\n')
	count_sites = 0
	while count_sites <= max_sites:
		snp = Snp(snpfiles, samples)
		if snp.data is None:
			break
		else:
			count_sites += 1
			count_samples = prevalence(snp, site_depth, freq=False)
			values = [snp.id()[0], snp.id()[1], snp.ref_allele(), count_samples]
			outfile.write('\t'.join([str(_) for _ in values])+'\n')

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
	batches = utility.batch_samples(samples, args['threads'])
	hq_snps = read_hq_snps(args['outdir'], args['species_id'])
	for index, batch in enumerate(batches):
		a = {}
		a['index'] = index
		a['samples'] = batch
		a['species_id'] = args['species_id']
		a['hq_snps'] = hq_snps
		a['outdir'] = tempdir
		list.append(a)
	utility.parallel(build_mini_matrix, list, args['threads'])
	# open temporary snp matrixes
	infiles = {}
	for type in ['ref_freq', 'depth', 'alt_allele']:
		files = []
		for index, batch in enumerate(batches):
			inpath = '%s/%s.snps.%s.%s' % (tempdir, args['species_id'], type, index)
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
	for type in ['ref_freq', 'depth', 'alt_allele']:
		for file in infiles[type]: file.close()
	shutil.rmtree(tempdir, ignore_errors=True)

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
	outfile = open('%s/%s.snps.fasta' % (args['outdir'], args['species_id']), 'w')
	hq_snps = read_hq_snps(args['outdir'], args['species_id']) # get snp list
	n = len(hq_snps)
	for sample in samples:
		i = 0 # index position in hq_snps
		outfile.write('>%s\n' % sample.id) # sequence header
		inpath = '%s/snps/snps/%s.snps.gz' % (sample.dir, args['species_id'])
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
	utility.print_copyright()
	print_arguments(args)

	print("Loading samples")
	samples = identify_samples(args)
	
	if args['snps']:
		#print("Identifying and writing hq snps") # TESTING
		#identify_snps(args, samples)
		print("Annotating snps")
		annotate_snps.main(args)

	if args['freq']:
		print("Writing allele frequencies & depths")
		build_snp_matrix(args, samples)

	if args['cons']:
		print("Writing consensus sequences")
		write_consensus(args, samples)

	if args['tree']:
		print("Building phylogenetic tree")
		build_tree(args)

			
