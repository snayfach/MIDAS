#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

## TO DO
# build trees
# add other references from genome-cluster
# write snp statistics (depth, prevalence, diversity, etc.)

__version__ = '0.0.2'

import argparse, sys, os, gzip, subprocess, tempfile

def parse_arguments():
	""" Parse command line arguments """
	
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
		
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', '--indir', type=str, dest='in', help='input directory', required=True)
	io.add_argument('-o', '--outdir', type=str, dest='out', help='output directory', required=True)
	io.add_argument('-g', '--genome_cluster', type=str,  help='genome cluster id', required=True)
	io.add_argument('-m', '--matrix', type=str,  help='reference SNP matrix')
	io.add_argument('-t', '--tree', default=False, action='store_true',  help='build phylogenetic tree')
	
	sample = parser.add_argument_group('Sample filters')
	sample.add_argument('--sample_list', dest='sample_list', type=str,
		default=None, help='file of sample ids to include; each line should contain one id')
	sample.add_argument('--sample_depth', dest='sample_depth', type=float,
		default=2.0, help='min average read depth per sample (2.0)')
	sample.add_argument('--ref_coverage', dest='ref_coverage', type=float,
		default=0.4, help='min coverage of reference genome per sample (0.4)')
				
	snps = parser.add_argument_group('SNP filters')
	snps.add_argument('--snp_prev', dest='min_prev', type=float,
		default=1.0, help='min fraction of samples that contain SNP (1.0)')
	snps.add_argument('--snp_depth', dest='snp_depth', type=int,
		default=1.0, help='min # of reads supporting SNP per sample (1)')
	snps.add_argument('--no_fixed', action='store_true', default=False,
		help='exclude SNPs with the same consensus allele across samples (False)')
	snps.add_argument('--max_snps', dest='max_snps', type=int,
		default=float('Inf'), help='only use <= MAX_SNPS (use all)')
					
	args = vars(parser.parse_args())
	args['consensus'] = '%s/%s.fasta' % (args['out'], args['genome_cluster'])
	
	return args

def parse_snps_summary(inpath):
	""" Read in summary snp statistics for genome-clusters """
	snps_summary = {}
	infile = open(inpath)
	fields = next(infile).rstrip().split()
	for line in open(inpath):
		values = line.rstrip().split()
		rec = dict([(i,j) for i,j in zip(fields, values)])
		snps_summary[rec['cluster_id']] = rec
	return snps_summary

def identify_samples(args):
	""" Identify samples with sufficient depth for target genome-cluster """
	samples = []
	sample_list = [x.rstrip() for x in open(args['sample_list']).readlines()] if args['sample_list'] else None
	for sample_id in os.listdir(args['in']):
		# read in summary stats for sample across genome-clusters
		indir = '/'.join([args['in'], sample_id])
		inpath = '/'.join([indir, 'snps_summary_stats.txt'])
		if not os.path.isfile(inpath):
			continue
		else:
			snps_summary = parse_snps_summary(inpath)
		# check whether sample passes QC
		if args['sample_list'] and sample_id not in sample_list:
			continue
		elif args['genome_cluster'] not in snps_summary:
			continue
		elif float(snps_summary[args['genome_cluster']]['average_depth']) < args['sample_depth']:
			print sample_id
			print snps_summary
			quit()
			continue
		elif float(snps_summary[args['genome_cluster']]['fraction_covered']) < args['ref_coverage']:
			continue
		elif not os.path.isfile('%s/snps/%s.snps.gz' % (indir, args['genome_cluster'])):
			continue
		# keep sample
		else:
			samples.append(sample_id)
	if len(samples) == 0:
		sys.exit("Error: no samples met selection criteria!")
	return samples

def parse_snps(inpath):
	""" Yields formatted records from SNPs output """
	infile = gzip.open(inpath)
	next(infile)
	fields = ['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
			  'count_alleles', 'count_ref', 'count_alt', 'depth', 'ref_freq']
	for line in infile:
		yield dict([(i,j) for i,j in zip(fields, line.rstrip().split())])

def open_infiles(args, samples):
	""" Open SNP files for genome-cluster across all samples """
	infiles = {}
	for sample_id in samples:
		inpath = '%s/%s/snps/%s.snps.gz' % (args['in'], sample_id, args['genome_cluster'])
		infiles[sample_id] = parse_snps(inpath)
	return infiles

def open_outfiles(args, samples):
	""" Open output files and write headers """
	# open outfiles
	outfiles = {}
	outfiles['cons_seqs'] = open('%s/%s.fasta' % (args['out'], args['genome_cluster']), 'w')
	outfiles['freq_matrix'] = open('%s/%s.ref_freq' % (args['out'], args['genome_cluster']), 'w')
	outfiles['depth_matrix'] = open('%s/%s.depth' % (args['out'], args['genome_cluster']), 'w')
	outfiles['snp_list'] = open('%s/%s.hq_snps' % (args['out'], args['genome_cluster']), 'w')
	# write headers
	header = ['ref_id', 'ref_pos'] + samples
	outfiles['freq_matrix'].write('\t'.join(header)+'\n')
	outfiles['depth_matrix'].write('\t'.join(header)+'\n')
	return outfiles

def fetch_snp(snpfiles, samples):
	""" Fetch SNP data across samples """
	snp = []
	for sample_id in samples:
		snp.append(next(snpfiles[sample_id]))
	return snp

def is_fixed(snp, pass_qc):
	""" Determine if SNP has the same consensus allele across all samples that passed QC """
	cons_alleles = [snp_i['cons_allele'] for is_hq, snp_i in zip(pass_qc, snp) if is_hq]
	if all([allele == cons_alleles[0] for allele in cons_alleles]):
		return True
	else:
		return False

def snp_qc(snp, args):
	""" Determine whether nor not to keep snp for each sample """
	i = []
	for s in snp:
		if float(s['depth']) < args['snp_depth']:
			i.append(False)
		else:
			i.append(True)
	return i

def prevalence(list):
	""" Count fraction of True in list """
	return float(sum(list))/len(list)

def write_records(snp, outfiles, args):
	""" Write snp to outfiles """
	id = [snp[0]['ref_id'], snp[0]['ref_pos']]
	outfiles['freq_matrix'].write('\t'.join(id+[str(s['ref_freq']) for s in snp])+'\n')
	outfiles['depth_matrix'].write('\t'.join(id+[str(s['depth']) for s in snp])+'\n')
	outfiles['snp_list'].write('\t'.join(id)+'\n')

def filter_snps_write_matrices(args, snpfiles, outfiles, samples):
	""" Identify SNPs that pass QC and optionally write data for these to output matrices """
	filtered_snps = []
	inpath = '%s/%s/snps/%s.snps.gz' % (args['in'], snpfiles.keys()[0], args['genome_cluster'])
	dummyfile = gzip.open(inpath)
	next(dummyfile)
	for line in dummyfile: # outer loop: stop when eof reached
		snp = fetch_snp(snpfiles, samples) # inner loop: get data for snp across all samples
		pass_qc = snp_qc(snp, args) # True/False depending if sample passed/failed QC for SNP
		if prevalence(pass_qc) < args['min_prev']: # skip snps that do not occur in enough samples
			continue
		elif args['no_fixed'] and is_fixed(snp, pass_qc): # skip snps that are are fixed across samples that passed QC
			continue
		else:
			filtered_snps.append([snp[0]['ref_id'], snp[0]['ref_pos']]) # keep track of snps that passed qc
			write_records(snp, outfiles, args) # write data to output files
		if len(filtered_snps) == args['max_snps']: # stop after reaching max_snps (for testing only)
			break
	for file in snpfiles.values(): # close open file handles
		file.close()
	for file in outfiles.values(): # close open file handles
		file.close()

def write_consensus(args, samples):
	""" Write consensus sequences from samples """
	outfile = open(args['consensus'], 'w')
	
	filtered_snps = read_filtered_snps(args) # get snp list
	n = len(filtered_snps) # number of snps in list
	
	# write consensus seqs for samples
	for sample_id in samples:
		i = 0 # index position in snp list
		outfile.write('>%s\n' % sample_id) # sequence header
		inpath = '%s/%s/snps/%s.snps.gz' % (args['in'], sample_id, args['genome_cluster'])
		for snp in parse_snps(inpath):
			if i == n: # no more snps in list
				break
			elif [snp['ref_id'], snp['ref_pos']] != filtered_snps[i]: # snp not in list
				continue
			else: # snp in list
				outfile.write('-' if snp['cons_allele'] == 'NA' else snp['cons_allele']) # write base
				i += 1
		outfile.write('\n')
#	# add consensus seqs for reference
#	i = 0 # index position in snp list
#	outfile.write('>representative_genome\n') # sequence header
#	inpath = '%s/%s/snps/%s.snps.gz' % (args['in'], samples[0], args['genome_cluster'])
#	for snp in parse_snps(inpath):
#		if i == n: # no more snps in list
#			break
#		elif [snp['ref_id'], snp['ref_pos']] != filtered_snps[i]: # snp not in list
#			continue
#		else: # snp in list
#			outfile.write(snp['ref_allele']) # write base
#			i += 1
#	outfile.write('\n')


def read_filtered_snps(args):
	snps = []
	inpath = '%s/%s.hq_snps' % (args['out'], args['genome_cluster'])
	for line in open(inpath):
		snps.append(line.rstrip().split())
	return snps

def add_ref_alleles(args):
	""" Add partial genotypes from reference genomes """
	
	meta_snps = read_filtered_snps(args) # get snp list
	n_snps = len(meta_snps) # length of snp list
	ref_snps = gzip.open(args['matrix']) # open file of reference genotypes
	max_genomes = 300 # maximum genomes to use from matrix
	
	# open tempfiles for reference genotypes
	tempfiles = [] #
	for line in ref_snps:
		genomes = line.rstrip().split()[2:max_genomes+1]
		for genome in genomes:
			file = open(tempfile.mkstemp()[1], 'w')
			file.write('>%s\n' % genome)
			tempfiles.append(file)
		break
		
	# write reference genotypes to temp files
	i = 0 # index position in snp list
	for line in ref_snps: # get partial genotypes for reference genomes
		ref_snp = line.rstrip().split()
		if i == n_snps: # no more snps in meta_snp list
			break
		elif ref_snp[0:2] != meta_snps[i]: # ref_snp not in meta_snp list
			continue
		else: # ref_snp in meta_snp list
			genotypes = ref_snp[2:max_genomes+1]
			for genotype, outfile in zip(genotypes, tempfiles):
				outfile.write(genotype)
			i += 1 # increment index
	
	# append reference genotypes to fasta
	outfile = open(args['consensus'], 'a')
	for f in tempfiles:
		f.close()
		f = open(f.name)
		outfile.write(f.read()+'\n')
		os.remove(f.name)

def build_tree(args):
	"""	Use FastTree to build phylogenetic tree of consensus sequences """
	inpath = '%s/%s.fasta' % (args['out'], args['genome_cluster'])
	outpath = '%s/%s.tree' % (args['out'], args['genome_cluster'])
	p = subprocess.Popen('FastTree -nt -boot 100 < %s > %s' % (inpath, outpath), shell=True)
	p.wait()

if __name__ == '__main__':

	args = parse_arguments()
	if not os.path.isdir(args['out']): os.mkdir(args['out'])

	print("Identifying samples")
	samples = identify_samples(args)
	
	print("Opening input files")
	snpfiles = open_infiles(args, samples)
	
	print("Opening output files")
	outfiles = open_outfiles(args, samples)

	print("Identifying and writing hq snps")
	filter_snps_write_matrices(args, snpfiles, outfiles, samples)

	print("Writing consensus sequences")
	write_consensus(args, samples)

	if args['matrix']:
		print("Adding references sequences")
		add_ref_alleles(args)

	if args['tree']:
		print("Building phylogenetic tree")
		build_tree(args)

			
