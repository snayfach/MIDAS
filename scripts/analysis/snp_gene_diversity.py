#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, gzip, numpy as np
from midas import utility, analyze_snps as analyze

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: nucleotide_diversity.py outdir [options]

Description: 
	use core-genome SNPs to estimate population genetic parameters of genes within and between samples
	within-sample parameters are estimated using the frequency of SNPs within individual samples
	between-sample parameters are estimated using the the mean frequency of SNPs across samples
Input: 
	SNP allele frequency matrix (SPECIES.snps.ref_freq), SNP annotation file (SPECIES.snps.info)
Output: 
	per-gene nucleotide diversity, within and between samples
	per-gene non-synonymous to synonymous SNP ratio (pN/pS), within and between samples
""",
		epilog="""Examples:
1) Run using defaults:
nucleotide_diversity.py /path/to/outdir --ref_freq /path/to/SPECIES.snps.ref.freq --snp_info /path/to/SPECIES.snps.info

2) Use alternate definition of a SNP as 1% minor allele frequency:
nucleotide_diversity.py /path/to/outdir --maf 0.01 --ref_freq /path/to/SPECIES.snps.ref.freq --snp_info /path/to/SPECIES.snps.info

3) run a quick test:
nucleotide_diversity.py /path/to/outdir --max_sites 10000 --ref_freq /path/to/SPECIES.snps.ref.freq --snp_info /path/to/SPECIES.snps.info

""")
	parser.add_argument('--out', type=str,
		help="""path to output file""")
	parser.add_argument('--inbase', metavar='PATH', type=str, required=True,
		help="""basename for input files""")
	parser.add_argument('--max_sites', metavar='INT', type=int, required=False,
		help="""maximum number of sites to use from input matrix
useful for quick tests (use all)""")
	parser.add_argument('--maf', type=float, metavar='FLOAT', required=False, default=0.05,
		help="""minor allele frequency for SNP definition (0.05)""")
	parser.add_argument('--samples', type=str, metavar='STR', required=False,
		help="""comma-separated list of samples to use for computing diversity metrics (use all)""")
	args = vars(parser.parse_args())
	args.update(st.init_paths(args))
	return args

def compute_diversity(args, samples):
	""" Compute number of sites and snps """
	# count sites and snps
	for site in st.parse_sites(args):
		freqs = []
		gene_id = site.info['gene_id']
		if gene_id == 'NA': continue # skip sites outside genes
		# within-sample
		for sample_id, freq, count in site.iterate(samples):
			freqs.append(freq)
			samples[sample_id].init_gene(gene_id) # initialize gene
			samples[sample_id].genes[gene_id].update_stats(site.info, freq, args['maf']) # record info
		# between-sample
		if len(freqs) > 0:
			samples['between'].init_gene(gene_id) # initialize gene
			samples['between'].genes[gene_id].update_stats(site.info, freq, args['maf']) # record info
	# compute diversity metrics
	for sample in samples.values():
		sample.compute_diversity()
		for gene in sample.genes.values():
			gene.compute_diversity()

def write_per_gene_results(args, samples):
	""" Write per-gene results to output file(s) """
	outfile = open(args['out'], 'w')
	fields = ['sample_id', 'gene_id', 'pi_bp', 'pnps_snps', 'pnps_pi', 'sites', 'snps',
			  'sites_NC', 'sites_1D', 'sites_2D', 'sites_3D', 'sites_4D',
			  'snps_NC', 'snps_1D', 'snps_2D', 'snps_3D', 'snps_4D',
			  'pi_NC', 'pi_1D', 'pi_2D', 'pi_3D', 'pi_4D']
	outfile.write('\t'.join(fields)+'\n')
	gene_ids = set([])
	for sample in samples.values():
		for gene_id in sample.genes:
			gene_ids.add(gene_id)
	for sample in samples.values():
		for gene_id in gene_ids:
			if gene_id not in sample.genes:
				gene = Gene(gene_id)
			else:
				gene = sample.genes[gene_id]
			record = gene.format_record(sample.id, fields)
			outfile.write(record)

if __name__ == '__main__':
	args = parse_arguments()
	utility.print_copyright()
	
	print("Loading samples...")
	samples = st.init_samples(args)
	
	print("Estimating diversity metrics...")
	compute_diversity(args, samples)

	print("Writing results...")
	write_per_gene_results(args, samples)

