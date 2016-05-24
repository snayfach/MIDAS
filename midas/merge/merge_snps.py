#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, shutil, numpy as np
from midas import utility
from midas.merge import merge
from midas import analyze_snps

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

def open_matrices(outbase, sample_ids, index=None):
	""" Open matrices and write headers """
	matrices = {}
	for type in ['ref_freq', 'depth', 'alt_allele']:
		if index is None: outpath = '%s.snps.%s' % (outbase, type)
		else: outpath = '%s.snps.%s.%s' % (outbase, type, index)
		matrices[type] = open(outpath, 'w')
		matrices[type].write('\t'.join(['site_id']+sample_ids)+'\n')
	return matrices

def build_snp_matrix(species_id, samples, args):
	""" Split up samples into batches, merge each batch, merge together batches """
	# build temp matrixes in parallel
	list = []
	tempdir = '%s/%s/temp' % (args['outdir'], species_id)
	if not os.path.isdir(tempdir): os.mkdir(tempdir)
	batches = utility.batch_samples(samples, args['threads'])
	for index, batch in enumerate(batches):
		list.append({'index':index, 'samples':batch, 'species_id':species_id,
		            'tempdir':tempdir, 'max_sites':args['max_sites']})
	utility.parallel(temp_matrix, list, args['threads'])
	# merge temp matrixes
	merge_matrices(tempdir, species_id, samples, batches, args)

def temp_matrix(tempdir, species_id, samples, index, max_sites):
	""" Build SNP matrices using a subset of total samples """
	sample_ids = [s.id for s in samples]
	outbase = '%s/%s' % (tempdir, species_id)
	matrices = open_matrices(outbase, sample_ids, index)
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
			site_id = '|'.join([records[0]['ref_id'], records[0]['ref_pos']])
			for field in ['ref_freq', 'depth', 'alt_allele']:
				values = [rec[field] for rec in records]
				matrices[field].write(site_id+'\t'+'\t'.join(values)+'\n')

def merge_matrices(tempdir, species_id, samples, batches, args):
	""" Merge together temp SNP matrices """
	if len(batches) == 1: # if only one batch, just rename files
		for type in ['ref_freq', 'depth', 'alt_allele']:
			inpath = '%s/%s.snps.%s.0' % (tempdir, species_id, type)
			outpath = '%s/%s.snps.%s' % (tempdir, species_id, type)
			shutil.move(inpath, outpath)
	else:  # if > one batch, merge temp matrices
		# open temporary matrixes
		infiles = {}
		for type in ['ref_freq', 'depth', 'alt_allele']:
			files = []
			for index, batch in enumerate(batches):
				inpath = '%s/%s.snps.%s.%s' % (tempdir, species_id, type, index)
				files.append(open(inpath))
			infiles[type] = files
		# merge temporary matrixes
		sample_ids = [s.id for s in samples]
		outbase = '%s.%s' % (tempdir, species_id)
		matrices = open_matrices(outbase, sample_ids)
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

#def build_tree(args):
#	"""	Use FastTree to build phylogenetic tree of consensus sequences """
#	inpath = '%s/%s.fasta' % (args['outdir'], args['species_id'])
#	outpath = '%s/%s.tree' % (args['outdir'], args['species_id'])
#	p = subprocess.Popen('FastTree -nt -boot 100 < %s > %s' % (inpath, outpath), shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
#	out, err = p.communicate()

def filter_snp_matrix(species_id, samples, args):
	""" Extract subset of site from SNP-matrix """
	tempbase = '%s/%s/temp/%s' % (args['outdir'], species_id, species_id)
	outbase = '%s/%s/%s' % (args['outdir'], species_id, species_id)
	if args['site_prev'] == 0:
		for type in ['ref_freq', 'depth', 'alt_allele']:
			shutil.move('%s.snps.%s' % (tempbase,type), '%s.snps.%s' % (outbase,type))
	else:
		outdir = os.path.join(args['outdir'], species_id)
		sample_ids = [s.id for s in samples]
		matrices = open_matrices(outbase, sample_ids)
		siteinfo = open('%s/%s.snps.info' % (outdir, species_id), 'w')
		siteinfo.write('\t'.join(['site_id', 'site_prev', 'mean_ref_freq', 'alt_props'])+'\n')
		for site in analyze_snps.parse_sites(tempbase, site_depth=args['site_depth'], max_sites=args['max_sites']):
			if (site.prev >= args['site_prev'] # site prevalence filter
					and (args['site_maf'] == 0.0 # site minor allele frequency filter
							or (site.maf >= args['site_maf']))):
				site.write_info(siteinfo, args)
				matrices['ref_freq'].write(site.id+'\t'+'\t'.join(site.freqs)+'\n')
				matrices['depth'].write(site.id+'\t'+'\t'.join(site.depths)+'\n')
				matrices['alt_allele'].write(site.id+'\t'+'\t'.join(site.alleles)+'\n')

def run_pipeline(args):

	print("Identifying species")
	species = merge.select_species(args, type='snps')
	for sp in species:

		print "Merging: %s (id:%s) for %s samples" % (sp.consensus_name, sp.id, len(sp.samples))

		print("  merging per-sample statistics")
		merge.write_summary_stats(sp.id, sp.samples, args, 'snps')

		print("  merging per-site statistics")
		build_snp_matrix(sp.id, sp.samples, args)

		print("  extracting specified sites")
		filter_snp_matrix(sp.id, sp.samples, args)

#		print("  removing temporary files")
#		shutil.rmtree('%s/%s/temp' % (args['outdir'], species_id))




			
