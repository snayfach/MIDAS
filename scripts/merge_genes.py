#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.2'

import argparse, sys, os, gzip
from collections import defaultdict

def print_copyright():
	print ("")
	print ("PhyloCNV: species abundance and strain-level genomic variation from metagenomes")
	print ("version %s; github.com/snayfach/PhyloCNV" % __version__)
	print ("Copyright (C) 2015 Stephen Nayfach")
	print ("Freely distributed under the GNU General Public License (GPLv3)")
	print ("")

def parse_arguments():
	""" Parse command line arguments """
	
	parser = argparse.ArgumentParser(
		usage='%s [options]' % os.path.basename(__file__),
		description="""Merge gene copy-number variants for an individual species across samples.
			Outputs include: a gene copy-number matrix, a gene presence/absence matrix, and a gene read-depth matrix""")
	parser.add_argument('-v', '--verbose', default=False, action='store_true', help='verbose')
	
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""input to results from 'run_phylo_cnv.py genes'.
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
			a map of species ids to species names can be found in 'ref_db/annotations.txt'""")
	io.add_argument('-o', type=str, dest='outdir', required=True,
		help="""output directory.
			output files: <outdir>/<species_id>.copynum, <outdir>/<species_id>.presabs, <outdir>/<species_id>.gene_depth""")
	
	sample = parser.add_argument_group('Sample filters\n(determine which samples are included in output)')
	sample.add_argument('--marker_coverage', type=float, default=3.0,
		help="""minimum read depth per sample across 15 phylogenetic marker genes (3.0)""")
	sample.add_argument('--gene_coverage', type=float, default=0.0,
		help="""minimum read depth per sample across all genes with non-zero coverage (0.0)""")
	sample.add_argument('--max_samples', type=int,
		help="""maximum number of samples to process. useful for testing (use all)""")

	gene = parser.add_argument_group('Presence/Absence')
	gene.add_argument('--min_copy', type=float, default=0.35,
		help="""genes >= MIN_COPY are classified as present and genes < MIN_COPY are classified as absent (0.35)""")
	gene.add_argument('--cluster_pid', type=str, dest='cluster_pid', default='95', choices=['75', '80', '85', '90', '95', '99'],
		help="""gene family percent identity. small values: fewer, larger gene families. large values: more, smaller gene families (95)""")
	gene.add_argument('--add_ref', default=False, action='store_true',
		help="""include gene presence/absence for reference genomes""")
	
	args = vars(parser.parse_args())
	args['db'] = '%s/ref_db/genome_clusters' % os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	
	return args

def print_arguments(args):
	print ("-------------------------------------------------------")
	print ("Merge Genes Parameters")
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Species identifier: %s" % args['species_id'])
	print ("Sample selection criteria:")
	if args['marker_coverage']:
		print ("  >=%s average coverage across 15 universal-single-copy genes" % args['marker_coverage'])
	if args['gene_coverage']:
		print ("  >=%s average coverage across all genes with non-zero coverage" % args['gene_coverage'])
	if args['max_samples']:
		print ("  analyze up to %s samples" % args['max_samples'])
	print ("Gene presence/absence criterea:")
	print ("  present (1): genes with copy number >=%s" % args['min_copy'])
	print ("  absent (0): genes with copy number <%s" % args['min_copy'])
	print ("  cluster genes at %s percent identity" % args['cluster_pid'])
	if args['add_ref']:
		print ("  adding gene presence/absences for reference genomes")
	print ("-------------------------------------------------------")
	print ("")

def read_file(inpath):
	""" Read in summary snp statistics for genome-clusters """
	d = {}
	infile = open(inpath)
	fields = next(infile).rstrip().split()
	for line in open(inpath):
		values = line.rstrip().split()
		rec = dict([(i,j) for i,j in zip(fields, values)])
		d[rec['cluster_id']] = rec
	return d

def list_samples(input, intype):
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

def identify_samples(args):
	""" Identify samples with sufficient depth for target genome-cluster """
	samples_dirs = []
	for sample_dir in list_samples(args['input'], args['intype']):
		# read in summary stats for sample across genome-clusters
		inpath = '/'.join([sample_dir, 'genes_summary_stats.txt'])
		sample_id = os.path.basename(sample_dir)
		if not os.path.isfile(inpath):
			print("  warning: no data for sample_id %s" % sample_id)
			continue
		genes_summary = read_file(inpath)
		# check whether sample passes QC
		if args['species_id'] not in genes_summary:
			continue
		elif float(genes_summary[args['species_id']]['phyeco_coverage']) < args['marker_coverage']:
			continue
		elif float(genes_summary[args['species_id']]['mean_coverage']) < args['gene_coverage']:
			continue
		# sample passes qc
		else:
			samples_dirs.append(sample_dir)
			# only keep max_samples if specified
			if args['max_samples'] and len(samples_dirs) >= args['max_samples']:
				break
	if len(samples_dirs) == 0:
		sys.exit("Error: species_id failed to pass quality-control in all samples")
	return samples_dirs

def parse_genes(inpath):
	""" Yields formatted records from coverage output """
	infile = gzip.open(inpath)
#	fields = next(infile).rstrip().split()
	next(infile)
	fields = ['gene_id', 'raw_coverage', 'normalized_coverage']
	for line in infile:
		yield dict([(i,j) for i,j in zip(fields, line.rstrip().split())])

def read_fam_map(args):
	""" Read gene family map """
	fam_map = {}
	inpath = '%s/%s/gene_family_map.txt.gz' % (args['db'], args['species_id'])
	infile = gzip.open(inpath)
	fields = next(infile).rstrip().split()
	for line in infile:
		values = line.rstrip().split()
		map = dict([(f,v) for f,v in zip(fields, values)])
		fam_map[map['99']] = map[args['cluster_pid']]
	return fam_map

def compute_sample_copy_num(sample_dirs, args, fam_map):
	""" Compute gene copy numbers for samples """
	fam_copy_num = {}
	for sample_dir in sample_dirs:
		sample_id = os.path.basename(sample_dir)
		copy_num = defaultdict(float)
		inpath = '%s/coverage/%s.cov.gz' % (sample_dir, args['species_id'])
		for r in parse_genes(inpath):
			fam_id = fam_map[r['gene_id']]
			copy_num[fam_id] += float(r['normalized_coverage'])
		fam_copy_num[sample_id] = copy_num
	return fam_copy_num

def compute_sample_read_depth(sample_dirs, args, fam_map):
	""" Compute gene read depth for samples """
	fam_read_depth = {}
	for sample_dir in sample_dirs:
		sample_id = os.path.basename(sample_dir)
		read_depth = defaultdict(float)
		inpath = '%s/coverage/%s.cov.gz' % (sample_dir, args['species_id'])
		for r in parse_genes(inpath):
			fam_id = fam_map[r['gene_id']]
			read_depth[fam_id] += float(r['raw_coverage'])
		fam_read_depth[sample_id] = read_depth
	return fam_read_depth

def read_genome_ids(args):
	""" Read in genome ids for genome cluster """
	genome_ids = set([])
	inpath = '%s/%s/genomes.txt.gz' % (args['db'], args['species_id'])
	infile = gzip.open(inpath); next(infile)
	for line in infile:
		genome_id = line.split('\t')[1]
		genome_ids.add(genome_id)
	return list(genome_ids)

def compute_ref_copy_num(genome_ids, args, fam_map):
	""" Compute gene copy numbers for reference genomes """
	ref_copy_num = dict([(_,defaultdict(float)) for _ in genome_ids])
	inpath = '%s/%s/ref_copy_num.txt.gz' % (args['db'], args['species_id'])
	infile = gzip.open(inpath); next(infile)
	for line in infile:
		genome_id, gene_id, copy_number = line.rstrip().split()
		fam_id = fam_map[gene_id]
		ref_copy_num[genome_id][fam_id] += float(copy_number)
	return ref_copy_num

def presence_absence(copy_nums, min_copynum):
	""" Compute presence/absences from normlized coverage values """
	presabs = [1 if copy_num >= min_copynum else 0 for copy_num in copy_nums]
	return presabs

def write_with_ref(sample_ids, genome_ids, sample_copy_num, ref_copy_num, args):
	# open outfiles
	outfiles = {}
	outfiles['copynum'] = open('%s/%s.copynum' % (args['outdir'], args['species_id']), 'w')
	outfiles['presabs'] = open('%s/%s.presabs' % (args['outdir'], args['species_id']), 'w')
	# write headers
	header = ['gene_id'] + sample_ids + genome_ids
	outfiles['copynum'].write('\t'.join(header)+'\n')
	outfiles['presabs'].write('\t'.join(header)+'\n')
	# write values
	for gene_id in genes:
		copy_nums = []
		outfiles['copynum'].write('%s\t' % gene_id) # write gene id
		outfiles['presabs'].write('%s\t' % gene_id)
		for sample_id in sample_ids:
			try: copy_nums.append(sample_copy_num[sample_id][gene_id])
			except: copy_nums.append(0.0)
		for genome_id in genome_ids:
			try: copy_nums.append(ref_copy_num[genome_id][gene_id])
			except: copy_nums.append(0.0)
		presabs = presence_absence(copy_nums, args['min_copy'])
		outfiles['copynum'].write('\t'.join([str(_) for _ in copy_nums])+'\n') # write gene values
		outfiles['presabs'].write('\t'.join([str(_) for _ in presabs])+'\n')

def write_wo_ref(sample_ids, sample_copy_num, args):
	# open outfiles
	outfiles = {}
	outfiles['copynum'] = open('%s/%s.copynum' % (args['outdir'], args['species_id']), 'w')
	outfiles['presabs'] = open('%s/%s.presabs' % (args['outdir'], args['species_id']), 'w')
	# write headers
	header = ['gene_id'] + sample_ids
	outfiles['copynum'].write('\t'.join(header)+'\n')
	outfiles['presabs'].write('\t'.join(header)+'\n')
	# write values
	for gene_id in genes:
		copy_nums = []
		outfiles['copynum'].write('%s\t' % gene_id) # write gene id
		outfiles['presabs'].write('%s\t' % gene_id)
		for sample_id in sample_ids:
			try: copy_nums.append(sample_copy_num[sample_id][gene_id])
			except: copy_nums.append(0.0)
		presabs = presence_absence(copy_nums, args['min_copy'])
		outfiles['copynum'].write('\t'.join([str(_) for _ in copy_nums])+'\n') # write gene values
		outfiles['presabs'].write('\t'.join([str(_) for _ in presabs])+'\n')

def write_read_depth(sample_ids, sample_read_depth, args):
	# open outfile
	outfile = open('%s/%s.gene_depth' % (args['outdir'], args['species_id']), 'w')
	outfile.write('\t'.join(['gene_id'] + sample_ids)+'\n')
	# write values
	for gene_id in genes:
		outfile.write(gene_id)
		for sample_id in sample_ids:
			read_depth = sample_read_depth[sample_id][gene_id]
			outfile.write('\t%s' % str(read_depth))
		outfile.write('\n')

if __name__ == '__main__':

	args = parse_arguments()
	if args['verbose']: print_copyright()
	if args['verbose']: print_arguments(args)
	if not os.path.isdir(args['outdir']): os.mkdir(args['outdir'])

	if args['verbose']: print("Identifying samples with species")
	sample_dirs = identify_samples(args) # id samples with sufficient depth
	sample_ids = [os.path.basename(_) for _ in sample_dirs]
	if args['verbose']: print("  %s samples with species" % len(sample_ids))
	
	if args['verbose']: print("Mapping gene ids")
	fam_map = read_fam_map(args) # map 99% gene ids to lower level
	genes = set(fam_map.values())
	
	if args['verbose']: print("Computing gene copy numbers for samples")
	sample_copy_num = compute_sample_copy_num(sample_dirs, args, fam_map) # gene copy numbers across samples
	sample_read_depth = compute_sample_read_depth(sample_dirs, args, fam_map)
	
	if args['add_ref']:
		if args['verbose']: print("Computing gene copy numbers for reference genomes")
		genome_ids = read_genome_ids(args)
		ref_copy_num = compute_ref_copy_num(genome_ids, args, fam_map)

	if args['verbose']: print("Writing results")
	write_read_depth(sample_ids, sample_read_depth, args)
	if args['add_ref']:
		write_with_ref(sample_ids, genome_ids, sample_copy_num, ref_copy_num, args)
	else:
		write_wo_ref(sample_ids, sample_copy_num, args)



