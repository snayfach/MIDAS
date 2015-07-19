#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3

__version__ = '0.0.1'

import argparse, sys, os
#phylo_cnv_path = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/phylo_cnv'
#sys.path.append(phylo_cnv_path)
from phylo_cnv import phylo_cnv

def parse_arguments():
	""" Parse command line arguments """
	
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
	
	parser.add_argument('--version', action='version', version='MicrobeCNV %s' % __version__)
	parser.add_argument('-v', '--verbose', action='store_true', default=False)
	
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-1', type=str, dest='m1', help='FASTQ file containing 1st mate if paired or unpaired reads')
	io.add_argument('-2', type=str, dest='m2', help='FASTQ file containing 2nd mate if paired')
	io.add_argument('-D', type=str, dest='db_dir', help='Directory of bt2 indexes for genome clusters')
	io.add_argument('-o', type=str, dest='out', help='Directory for output files')

	pipe = parser.add_argument_group('Pipeline')
	pipe.add_argument('--all', action='store_true', dest='all',
		default=False, help='Run entire pipeline')
	pipe.add_argument('--profile', action='store_true', dest='profile',
		default=False, help='Estimate genome-cluster abundance using MicrobeSpecies')
	pipe.add_argument('--align', action='store_true', dest='align',
		default=False, help='Align reads to genome-clusters')
	pipe.add_argument('--map', action='store_true', dest='map',
		default=False, help='Assign reads to mapping locations')
	pipe.add_argument('--cov', action='store_true', dest='cov',
		default=False, help='Compute coverage of pangenomes')
	pipe.add_argument('--extract', action='store_true', dest='extract',
		default=False, help='Extract mapped reads from bam file & write to FASTQ')
	pipe.add_argument('--remap', action='store_true', dest='remap',
		default=False, help='Re-map reads to representative genomes')
	pipe.add_argument('--snps', action='store_true', dest='snps',
		default=False, help='Run samtools mpileup & estimate SNP frequencies')
				
	gc = parser.add_argument_group('Genome-cluster inclusion (choose one)')
	gc.add_argument('--gc_topn', type=int, dest='gc_topn', default=5, help='Top N most abundant (5)')
	gc.add_argument('--gc_cov', type=float, dest='gc_cov', help='Coverage threshold (None)')
	gc.add_argument('--gc_rbun', type=float, dest='gc_rbun', help='Relative abundance threshold (None)')
	gc.add_argument('--gc_id', type=str, dest='gc_id', help='Identifier of specific genome cluster (None)')
	gc.add_argument('--gc_list', type=str, dest='gc_list', help='Comma-separated list of genome cluster ids (None)')
	
	reads = parser.add_argument_group('Read selection')
	reads.add_argument('--rd_ms', type=int, dest='reads_ms',
		default=5000000, help='# reads to use for estimating genome-cluster abundance (5,000,000)')
	reads.add_argument('--rd_align', type=int, dest='reads_align',
		help='# reads to use for pangenome alignment (All)')
	reads.add_argument('--rd_batch', type=int, dest='rd_batch', default=5000000,
		help='Batch size in # reads. Smaller batch sizes requires less memory, but can take longer to run (5,000,000)')

	map = parser.add_argument_group('Mapping')
	map.add_argument('--map_pid', type=float, dest='pid',
		default=93, help='Minimum percent identity between read and reference (93.0)')
		
	snps = parser.add_argument_group('SNP detection')
	snps.add_argument('--snps_mapq', type=str, dest='snps_mapq',
		default='0', help='Minimum map quality (0)')
	snps.add_argument('--snps_baseq', type=str, dest='snps_baseq',
		default='0', help='Minimum base quality (0)')
				
	mask = parser.add_argument_group('Leave-One-Out (for simulated data only)')
	mask.add_argument('--tax_mask', action='store_true', dest='tax_mask',
		default=False, help='Discard alignments for reads and ref seqs from the same genome')
	mask.add_argument('--tax_map', type=str, dest='tax_map',
		help='File mapping read ids to genome ids')
	
	return parser.parse_args()

if __name__ == '__main__':
	args = vars(parse_arguments())
	if not os.path.isdir(args['out']): os.mkdir(args['out'])
	phylo_cnv.run_pipeline(args)

