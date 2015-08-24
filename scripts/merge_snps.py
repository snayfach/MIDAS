#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

## TO DO
# build trees
# add other references from genome-cluster
# write snp statistics (depth, prevalence, diversity, etc.)

__version__ = '0.0.2'

import argparse, sys, os, gzip

def parse_arguments():
	""" Parse command line arguments """
	
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
		
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', '--indir', type=str, dest='in', help='Input directory', required=True)
	io.add_argument('-o', '--outdir', type=str, dest='out', help='Output directory', required=True)
	io.add_argument('-g', '--genome_cluster', type=str,  help='Genome cluster identifer', required=True)
	
	sample = parser.add_argument_group('Sample filters')
	sample.add_argument('--sample_depth', dest='sample_depth', type=int,
		default=2.0, help='Min mean read depth at sites with non-zero coverage (2)')
	sample.add_argument('--ref_coverage', dest='ref_coverage', type=float,
		default=0.4, help='Min fraction of reference sites with non-zero coverage (0.4)')
				
	snps = parser.add_argument_group('SNP filters')
	snps.add_argument('--snp_prev', dest='min_prev', type=float,
		default=1.0, help='Minimum fraction of samples that contain SNP (1.0)')
	snps.add_argument('--snp_depth', dest='snp_depth', type=int,
		default=1.0, help='Minimum number of reads per sample supporting SNP (1)')
	snps.add_argument('--no_fixed', action='store_true', default=False,
		help='Exclude SNPs with the same consensus allele across samples (False)')
	snps.add_argument('--max_snps', dest='max_snps', type=int,
		default=float('Inf'), help='Only use up to MAX_SNPS. Useful for testing (use all)')
		
	other = parser.add_argument_group('Other')
	other.add_argument('-v', '--verbose', action='store_true', default=False)
	other.add_argument('-t', '--threads', dest='threads', default=1, help='Number of threads to use')
			
	args = vars(parser.parse_args())
	
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
	for sample_id in os.listdir(args['in']):
		# read in summary stats for sample across genome-clusters
		indir = '/'.join([args['in'], sample_id])
		inpath = '/'.join([indir, 'snps_summary_stats.txt'])
		if not os.path.isfile(inpath):
			continue
		else:
			snps_summary = parse_snps_summary(inpath)
		# check whether sample passes QC
		if args['genome_cluster'] not in snps_summary:
			continue
		elif float(snps_summary[args['genome_cluster']]['average_depth']) < args['sample_depth']:
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
	return filtered_snps

def write_consensus(args, samples, filtered_snps):
	""" Write consensus sequences from samples and reference genome """
	outfile = open('%s/%s.fasta' % (args['out'], args['genome_cluster']), 'w')
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
	# add consensus seqs for reference
	i = 0 # index position in snp list
	outfile.write('>representative_genome\n') # sequence header
	inpath = '%s/%s/snps/%s.snps.gz' % (args['in'], samples[0], args['genome_cluster'])
	for snp in parse_snps(inpath):
		if i == n: # no more snps in list
			break
		elif [snp['ref_id'], snp['ref_pos']] != filtered_snps[i]: # snp not in list
			continue
		else: # snp in list
			outfile.write(snp['ref_allele']) # write base
			i += 1
	outfile.write('\n')

if __name__ == '__main__':

	args = parse_arguments()
	if not os.path.isdir(args['out']): os.mkdir(args['out'])

	# id samples with sufficient depth
	print("Identifying samples")
	samples = identify_samples(args)
	
	print("Opening input files")
	snpfiles = open_infiles(args, samples)
	
	print("Opening output files")
	outfiles = open_outfiles(args, samples)

	print("Identifying and writing hq snps")
	filtered_snps = filter_snps_write_matrices(args, snpfiles, outfiles, samples)

	print("Writing consensus sequences")
	write_consensus(args, samples, filtered_snps)

			
