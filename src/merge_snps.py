#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

## LIBRARIES
import sys, os, gzip, subprocess, tempfile, shutil
import annotate_snps, merge, utility

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
def write_samples(args, samples):
	""" Write sample list to file """
	outfile = open('%s/%s.snps.sample_ids' % (args['outdir'], args['species_id']), 'w')
	outfile.write('\n'.join([_.id for _ in samples])+'\n')
	outfile.close()
	
	outfile = open('%s/%s.snps.summary' % (args['outdir'], args['species_id']), 'w')
	outfile.write('\t'.join(['sample_id', 'average_depth', 'fraction_covered'])+'\n')
	for sample in samples:
		outfile.write(sample.id+'\t')
		outfile.write(str(sample.average_depth)+'\t')
		outfile.write(str(sample.fraction_covered)+'\n')
	outfile.close()

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
			sample.average_depth = stats[args['species_id']]['average_depth']
			sample.fraction_covered = stats[args['species_id']]['fraction_covered']
			samples.append(sample)
		if args['max_samples'] and len(samples) >= args['max_samples']:
			break
	if len(samples) == 0:
		sys.exit("\nError: no samples with sufficient coverage for species_id.\nTry running with more lenient parameters")
	else:
		write_samples(args, samples)
		print("  %s samples with species" % len(samples))
		return samples

def open_infiles(species_id, samples):
	""" Open SNP files for species across samples """
	infiles = {}
	for sample in samples:
		inpath = '%s/snps/snps/%s.snps.gz' % (sample.dir, species_id)
		infiles[sample.id] = utility.parse_file(inpath)
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
		if index is None: outpath = '%s/%s.snps.%s' % (outdir, species_id, type)
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
		for rec in utility.parse_file(inpath):
			if i == n: # no more snps in list
				break
			elif [rec['ref_id'], rec['ref_pos']] != hq_snps[i]: # snp not in list
				continue
			else: # snp in list
				outfile.write('-' if rec['cons_allele'] == 'NA' else rec['cons_allele']) # write base
				i += 1
		outfile.write('\n')

def build_tree(args):
	"""	Use FastTree to build phylogenetic tree of consensus sequences """
	inpath = '%s/%s.fasta' % (args['outdir'], args['species_id'])
	outpath = '%s/%s.tree' % (args['outdir'], args['species_id'])
	p = subprocess.Popen('FastTree -nt -boot 100 < %s > %s' % (inpath, outpath), shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	out, err = p.communicate()

def run_pipeline(args):

	print("Loading samples")
	samples = identify_samples(args)
	
	if args['snps']:
		print("Identifying and writing hq snps")
		identify_snps(args, samples)
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

			
