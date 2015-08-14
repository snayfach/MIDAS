#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3

__version__ = '0.0.1'

import argparse, sys, os
from phylo_cnv import phylo_cnv

def parse_arguments():
	""" Parse command line arguments """
	
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
	
	parser.add_argument('--version', action='version', version='MicrobeCNV %s' % __version__)
	parser.add_argument('-v', '--verbose', action='store_true', default=False)
	parser.add_argument('-t', '--threads', dest='threads', default=1, help='Number of threads to use')
	parser.add_argument('--tax_mask', dest='tax_mask', default=None, help=argparse.SUPPRESS)
	
	io = parser.add_argument_group('Input/Output (required)')
	io.add_argument('-1', type=str, dest='m1', help='FASTQ file containing 1st mate if paired or unpaired reads', required=True)
	io.add_argument('-2', type=str, dest='m2', help='FASTQ file containing 2nd mate if paired')
	io.add_argument('-D', type=str, dest='db_dir', help='Directory of bt2 indexes for genome clusters', required=True)
	io.add_argument('-o', type=str, dest='out', help='Directory for output files', required=True)

	pipe = parser.add_argument_group('Pipeline')
	pipe.add_argument('--all', action='store_true', dest='all',
		default=False, help='Run entire pipeline')
	pipe.add_argument('--species_profile', action='store_true', dest='species_profile',
		default=False, help='Estimate genome-cluster abundance')
	pipe.add_argument('--pangenome_build_db', action='store_true', dest='pangenome_build_db',
		default=False, help='Build bowtie2 database of pangenome centroids')
	pipe.add_argument('--pangenome_align', action='store_true', dest='pangenome_align',
		default=False, help='Align reads to genome-clusters')
	pipe.add_argument('--pangenome_cov', action='store_true', dest='pangenome_cov',
		default=False, help='Compute coverage of pangenomes')
	pipe.add_argument('--snps_build_db', action='store_true', dest='snps_build_db',
		default=False, help='Build bowtie2 database of representative genomes')
	pipe.add_argument('--snps_align', action='store_true', dest='snps_align',
		default=False, help='Align reads to representative genomes')
	pipe.add_argument('--snps_call', action='store_true', dest='snps_call',
		default=False, help='Run samtools mpileup & estimate SNP frequencies')

	profile = parser.add_argument_group('GC Abundance')
	profile.add_argument('--reads_gc', type=int, dest='reads_ms',
		default=5000000, help='# reads to use for estimating genome-cluster abundance (5000000)')

	gc = parser.add_argument_group('GC inclusion (choose one)')
	gc.add_argument('--gc_topn', type=int, dest='gc_topn', help='Top N most abundant (None)')
	gc.add_argument('--gc_cov', type=float, dest='gc_cov', help='Coverage threshold (None)')
	gc.add_argument('--gc_rbun', type=float, dest='gc_rbun', help='Relative abundance threshold (None)')
	gc.add_argument('--gc_id', type=str, dest='gc_id', help='Identifier of specific genome cluster or comma-separated list of ids (None)')
	
	map = parser.add_argument_group('Read Alignment/Mapping')
	map.add_argument('--align_speed', dest='align_speed',
		choices=[
			'very-fast', 'fast', 'sensitive', 'very-sensitive',
			'very-fast-local', 'fast-local', 'sensitive-local', 'very-sensitive-local'],
		help='alignment speed/sensitivity (sensitive)', default='sensitive')
	map.add_argument('--reads_align', type=int, dest='reads_align',
		help='# reads for pangenome or genome alignment (use all)')
	map.add_argument('--map_pid', type=float, dest='pid',
		default=93, help='Minimum percent ID between read and reference (93.0)')
	map.add_argument('--aln_cov', type=float, dest='aln_cov',
		default=0.70, help='Minimum alignment coverage of read (0.70)')
		
	snps = parser.add_argument_group('SNP detection')
	snps.add_argument('--snps_mapq', type=str, dest='snps_mapq',
		default='20', help='Minimum map quality (20)')
	snps.add_argument('--snps_baseq', type=str, dest='snps_baseq',
		default='20', help='Minimum base quality (20)')
				
	args = parser.parse_args()
	if args.tax_mask: args.tax_mask = args.tax_mask.split(',')
	if args.gc_id: args.gc_id = args.gc_id.split(',')
	
	return args

if __name__ == '__main__':
	args = vars(parse_arguments())
	if not os.path.isdir(args['out']): os.mkdir(args['out'])
	phylo_cnv.run_pipeline(args)

