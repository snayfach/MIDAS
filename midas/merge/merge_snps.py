#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, shutil, numpy as np
from midas import utility
from midas.merge import merge, annotate, snp_matrix

def store_data(snpfiles):
	""" List of records from specified sample_ids """
	x = []
	for file in snpfiles:
		try: x.append(next(file))
		except StopIteration: return None
	return x

def open_infiles(species_id, samples):
	""" Open SNP files for species across samples """
	infiles = []
	for sample in samples:
		inpath = '%s/snps/output/%s.snps.gz' % (sample.dir, species_id)
		infiles.append(utility.parse_file(inpath))
	return infiles

def open_matrices(outdir, sample_ids, index=None):
	""" Open matrices and write headers """
	matrices = {}
	for type in ['ref_freq', 'depth', 'alt_allele']:
		if index is None: outpath = '%s/snps_%s.txt' % (outdir, type)
		else: outpath = '%s/snps_%s.%s.txt' % (outdir, type, index)
		matrices[type] = open(outpath, 'w')
		matrices[type].write('\t'.join(['site_id']+sample_ids)+'\n')
	return matrices

def build_snp_matrix(species_id, samples, args):
	""" Split up samples into batches, merge each batch, merge together batches """
	# build temp matrixes in parallel
	list = []
	tempdir = '%s/%s/temp' % (args['outdir'], species_id)
	if not os.path.isdir(tempdir): os.mkdir(tempdir)
	batches = utility.batch_samples(samples, threads=1)
	for index, batch in enumerate(batches):
		temp_matrix(tempdir, species_id, batch, index, args['max_sites'])
	# merge temp matrixes
	merge_matrices(tempdir, species_id, samples, batches, args)

def temp_matrix(tempdir, species_id, samples, index, max_sites):
	""" Build SNP matrices using a subset of total samples """
	sample_ids = [s.id for s in samples]
	matrices = open_matrices(tempdir, sample_ids, index)
	snpfiles = open_infiles(species_id, samples)
	nsites = 0
	while True:
		records = store_data(snpfiles)
		if records is None: # eof
			break
		elif nsites >= max_sites:
			break
		else:
			nsites += 1
			site_id = '|'.join([records[0]['ref_id'], records[0]['ref_pos'], records[0]['ref_allele']])
			for field in ['ref_freq', 'depth', 'alt_allele']:
				values = [rec[field] for rec in records]
				matrices[field].write(site_id+'\t'+'\t'.join(values)+'\n')

def merge_matrices(tempdir, species_id, samples, batches, args):
	""" Merge together temp SNP matrices """
	if len(batches) == 1: # if only one batch, just rename files
		for type in ['ref_freq', 'depth', 'alt_allele']:
			inpath = '%s/snps_%s.0.txt' % (tempdir, type)
			outpath = '%s/snps_%s.txt' % (tempdir, type)
			shutil.move(inpath, outpath)
	else:  # if > one batch, merge temp matrices
		# open temporary matrixes
		infiles = {}
		for type in ['ref_freq', 'depth', 'alt_allele']:
			files = []
			for index, batch in enumerate(batches):
				inpath = '%s/snps_%s.%s.txt' % (tempdir, type, index)
				files.append(open(inpath))
			infiles[type] = files
		# merge temporary matrixes
		matrices = open_matrices(tempdir, sample_ids=[s.id for s in samples])
		for type in ['ref_freq', 'depth', 'alt_allele']:
			files = infiles[type]
			for file in files: next(file) # skip header
			while True:
				values = []
				try:
					for index, file in enumerate(files):
						v = next(file).rstrip().split()
						if index == 0: values += v
						else: values += v[1:]
					matrices[type].write('\t'.join(values)+'\n')
				except StopIteration:
					break
		for file in matrices.values(): file.close()
		# close temp files
		for type in ['ref_freq', 'depth', 'alt_allele']:
			for file in infiles[type]:
				file.close()

def format_dict(d):
	""" Format dictionary. ex: 'A:SYN|C:NS|T:NS|G:NS' """
	return '|'.join(['%s:%s' % (x, y) for x, y in d.items()])

def write_site_info(siteinfo, site_depth=None, site=None, header=None):
	""" Write site info to file """
	if header:
		fields = ['site_id', 'mean_freq', 'mean_depth', 'site_prev', 'allele_props', 'site_type', 'gene_id', 'amino_acids', 'snps']
		siteinfo.write('\t'.join(fields)+'\n')
	else:
		rec = []
		rec.append(site.id)
		rec.append(site.mean_freq())
		rec.append(site.mean_depth())
		rec.append(site.prev(site_depth))
		rec.append(format_dict(site.allele_props()))
		rec.append(site.site_type)
		rec.append(site.gene_id)
		rec.append(format_dict(site.amino_acids))
		rec.append(format_dict(site.snp_types))
		siteinfo.write('\t'.join([str(_) for _ in rec])+'\n')

def write_matrices(site, matrices):
	""" Write site to output matrices """
	matrices['ref_freq'].write(site.id+'\t'+'\t'.join(site.ref_freq)+'\n')
	matrices['depth'].write(site.id+'\t'+'\t'.join(site.depth)+'\n')
	matrices['alt_allele'].write(site.id+'\t'+'\t'.join(site.alt_allele)+'\n')

def filter_snp_matrix(species_id, samples, args):
	""" Extract subset of site from SNP-matrix """
	
	# init variables for site annotation
	gene_index = [0]
	contigs = annotate.read_genome(args['db'], species_id)
	genes = annotate.read_genes(args['db'], species_id, contigs)

	# open site matrixes
	outdir = os.path.join(args['outdir'], species_id)
	sample_ids = [s.id for s in samples]
	matrices = open_matrices(outdir, sample_ids)
	
	# open site info file & write header
	siteinfo = open('%s/snps_info.txt' % outdir, 'w')
	write_site_info(siteinfo, header=True)
	
	# parse genomic sites
	tempdir = '%s/%s/temp' % (args['outdir'], species_id)
	for index, site in enumerate(snp_matrix.parse_sites(tempdir)):
		if args['max_sites'] is not None and index >= args['max_sites']:
			break
		elif site.filter(args['site_depth'], args['site_prev'], args['site_maf']):
			continue
		else:
			annotate.annotate_site(site, genes, gene_index, contigs)
			write_site_info(siteinfo, args['site_depth'], site)
			write_matrices(site, matrices)

def write_readme(args, sp):
	outfile = open('%s/%s/readme.txt' % (args['outdir'], sp.id), 'w')
	outfile.write("""
Description of output files and file formats from 'merge_midas.py snps'

Output files
############
snps_ref_freq.txt  
  frequency of reference allele per genomic site and per sample (0.0)
snps_alt_allele.txt  
  alternate allele per genomic site and per sample
snps_depth.txt  
  number of reads mapped to genomic site per sample
snps_info.txt  
  metadata for genomic site
snps_summary.txt
  alignment summary statistics per sample
snps_log.txt
  log file containing parameters used

Output formats
############
snps_ref_freq.txt, snps_alt_allele.txt, snps_depth.txt,
  tab-delimited matrix files
  field names are sample ids
  row names are genome site ids
snps_summary.txt
  sample_id: sample identifier
  genome_length: number of base pairs in representative genome
  covered_bases: number of reference sites with at least 1 mapped read
  fraction_covered: proportion of reference sites with at least 1 mapped read
  mean_coverage: average read-depth across reference sites with at least 1 mapped read
snps_info.txt
  site_id: genomic site_id; format: ref_id|ref_pos|ref_allele
  mean_freq: average frequency of reference allele across samples
  mean_depth: average read-depth across samples
  site_prev: proportion of samples where site_id was covered with sufficient depth
  allele_props: pooled frequency of 4 nucleotides
  site_type: NC (non-coding), 1D, 2D, 3D, 4D (degeneracy)
  gene_id: gene that intersects site
  amino_acids: protein affect of all 4 possible nucleotides
  snps: SYN/NS for all 4 possible nucleotides

Additional information for species can be found in the reference database:
 %s/rep_genomes/%s
""" % (args['db'], sp.id) )
	outfile.close()

def merge_snps(args, species):
	log = open('%s/%s/snps_log.txt' % (args['outdir'], species.id), 'w')
	log.write("Merging: %s for %s samples\n" % (species.id, len(species.samples)))
	log.write("  merging per-sample statistics\n")
	merge.write_summary_stats(species.id, species.samples, args, 'snps')
	log.write("  merging per-site statistics\n")
	build_snp_matrix(species.id, species.samples, args)
	log.write("  extracting and annotating specified sites\n")
	filter_snp_matrix(species.id, species.samples, args)
	log.write("  removing temporary files\n")
	shutil.rmtree('%s/%s/temp' % (args['outdir'], species.id))
	log.close()
	write_readme(args, species)

def run_pipeline(args):

	print("Identifying species")
	species = merge.select_species(args, type='snps')
	
	print("Merging snps")
	batches =[]
	for species in species:
		batches.append({'args':args, 'species':species})
	utility.parallel(merge_snps, batches, args['threads'])





			
