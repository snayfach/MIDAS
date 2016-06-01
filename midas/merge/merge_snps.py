#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, shutil, numpy as np
from midas import utility
from midas.merge import merge, annotate_sites as annotate
from midas import analyze_snps as analyze

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
			site_id = '|'.join([records[0]['ref_id'], records[0]['ref_pos']])
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

#def build_tree(args):
#	"""	Use FastTree to build phylogenetic tree of consensus sequences """
#	inpath = '%s/%s.fasta' % (args['outdir'], args['species_id'])
#	outpath = '%s/%s.tree' % (args['outdir'], args['species_id'])
#	p = subprocess.Popen('FastTree -nt -boot 100 < %s > %s' % (inpath, outpath), shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
#	out, err = p.communicate()

def filter_site(site, args):
	""" Filter genome site based on prevalence and minor allele frequency """
	if site.prev < args['site_prev']:
		return True
	elif site.maf < args['site_maf']:
		return True
	else:
		return False

def format_dict(d):
	""" Format dictionary. ex: 'A:SYN|C:NS|T:NS|G:NS' """
	return '|'.join(['%s:%s' % (x, y) for x, y in d.items()])

def write_site_info(siteinfo, site=None, header=None):
	""" Write site info to file """
	if header:
		fields = ['site_id', 'mean_freq', 'mean_depth', 'site_prev', 'ref_allele',
			      'alt_alleles', 'site_type', 'gene_id', 'amino_acids', 'snps']
		siteinfo.write('\t'.join(fields)+'\n')
	else:
		rec = []
		rec.append(site.id)
		rec.append(site.mean_freq)
		rec.append(site.mean_depth)
		rec.append(site.prev)
		rec.append(site.ref_allele)
		rec.append(format_dict(site.alt_alleles()))
		rec.append(site.site_type)
		rec.append(site.gene_id)
		rec.append(format_dict(site.amino_acids))
		rec.append(format_dict(site.snp_types))
		siteinfo.write('\t'.join([str(_) for _ in rec])+'\n')

def write_matrices(site, matrices):
	""" Write site to output matrices """
	matrices['ref_freq'].write(site.id+'\t'+'\t'.join(site.freqs)+'\n')
	matrices['depth'].write(site.id+'\t'+'\t'.join(site.depths)+'\n')
	matrices['alt_allele'].write(site.id+'\t'+'\t'.join(site.alleles)+'\n')

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
	for site in analyze.parse_sites(tempdir, site_depth=args['site_depth'], max_sites=args['max_sites']):
		if filter_site(site, args):
			continue
		else:
			annotate.annotate_site(site, genes, gene_index, contigs)
			write_site_info(siteinfo, site)
			write_matrices(site, matrices)

def run_pipeline(args):

	print("Identifying species")
	species = merge.select_species(args, type='snps')
	for sp in species:

		print "Merging: %s (id:%s) for %s samples" % (sp.consensus_name, sp.id, len(sp.samples))

		print("  merging per-sample statistics")
		merge.write_summary_stats(sp.id, sp.samples, args, 'snps')

		print("  merging per-site statistics")
		build_snp_matrix(sp.id, sp.samples, args)

		print("  extracting and annotating specified sites")
		filter_snp_matrix(sp.id, sp.samples, args)

		print("  removing temporary files")
		shutil.rmtree('%s/%s/temp' % (args['outdir'], sp.id))




			
