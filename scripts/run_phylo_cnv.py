#!/usr/bin/env python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, platform
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from phylo_cnv import utility

def get_program():
	""" Get program specified by user (species, genes, or snps) """
	if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
		print('')
		print('Usage: run_phylo_cnv.py <command> [options]')
		print('')
		print('Note: use run_phylo_cnv.py <command> -h to view usage for a specific command')
		print('')
		print('Commands:')
		print('\tspecies\t estimate the abundance of 5,952 bacterial species')
		print('\tgenes\t identify gene copy number variants in abundant species')
		print('\tsnps\t identify single nucleotide variants in abundant species')
		quit()
	elif sys.argv[1] not in ['species', 'genes', 'snps']:
		sys.exit("Unrecognized command: '%s'" % sys.argv[1])
		quit()
	else:
		return sys.argv[1]

def get_arguments(program):
	""" Get arguments for specified program """
	if program == 'species':
		args = species_arguments()
	elif program == 'genes':
		args = gene_arguments()
	elif program == 'snps':
		args = snp_arguments()
	else:
		sys.error("Unrecognized program: '%s'" % program)
	args = utility.add_ref_db(args)
	return args

def check_arguments(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		check_species(args)
	elif program == 'genes':
		check_genes(args)
	elif program == 'snps':
		check_snps(args)
	else:
		sys.error("Unrecognized program: '%s'" % program)
	if platform.system() not in ['Linux', 'Darwin']:
		sys.exit("Operating system '%s' not supported" % system())

def print_arguments(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		print_species_arguments(args)
	elif program == 'genes':
		print_gene_arguments(args)
	elif program == 'snps':
		print_snp_arguments(args)
	else:
		sys.error("Unrecognized program: '%s'" % program)
	
def run_program(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		from phylo_cnv import species
		species.estimate_abundance(args)
	elif program == 'genes':
		from phylo_cnv import genes
		genes.run_pipeline(args)
	elif program == 'snps':
		from phylo_cnv import snps
		snps.run_pipeline(args)
	else:
		sys.error("Unrecognized program: '%s'" % program)


def species_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: run_phylo_cnv.py species outdir [options]

Description:
This script will map metagenomic reads to a database of phylogenetic marker genes using HS-BLASTN
Mapped reads are used estimate the read depth and relative abundance of bacterial species
Reads are mapped according to gene-specific, species-level mapping thresholds (94.5-98% DNA identity)
Reads that map equally well to 2 or more species are probabalistically assigned
""",
		epilog="""Examples:
1) run with defaults using a paired-end metagenome:
run_phylo_cnv.py species outdir -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz

2) run using a single-end metagenome with 4 CPUs and only 4M reads:
run_phylo_cnv.py species outdir -1 /path/to/reads_1.fq.gz -t 4 -n 4000000

3) run with exactly 80 base-pair reads:
run_phylo_cnv.py species outdir -1 /path/to/reads_1.fq.gz --read_length 80

4) quantify species abundance using a 16S database:
run_phylo_cnv.py species outdir -1 /path/to/reads_1.fq.gz --db_type ssuRNA
	""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str, help='Path to directory to store results. Name should correspond to sample identifier.')
	parser.add_argument('-1', type=str, dest='m1', required=True,
		help="FASTA/FASTQ file containing 1st mate if paired or unpaired reads")
	parser.add_argument('-2', type=str, dest='m2',
		help="FASTA/FASTQ file containing 2nd mate if paired")
	parser.add_argument('-n', type=int, dest='max_reads',
		help="""Number of reads to use from input file(s) (use all)""")
	parser.add_argument('-t', dest='threads', default=1,
		help="""Number of threads to use for database search (1)""")
	parser.add_argument('--db_type', choices=['phyeco', 'ssuRNA'], metavar='', default='phyeco',
		help="""Reference database. Choices:\n'phyeco': universal-single-copy protein family database (default)\n'ssuRNA': 16S ribosomal rna database""")
	parser.add_argument('-D', type=str, dest='db', help=argparse.SUPPRESS)
	parser.add_argument('--remove_temp', dest='remove_temp', default=False, action='store_true',
		help="Remove temporary files, including BLAST output")
	parser.add_argument('--word_size', type=int, metavar='INT', default=28,
		help="""Word size for BLAST search (28)\nUse word sizes > 16 for greatest efficiency.""")
	parser.add_argument('--mapid', type=float, metavar='FLOAT',
		help="""Discard reads with alignment identity < MAPID\nBy default gene-specific species-level cutoffs are used\nValues between 0-100 accepted""")
	parser.add_argument('--aln_cov', type=float, metavar='FLOAT', default=0.75,
		help="""Discard reads with alignment coverage < ALN_COV (0.75)\nValues between 0-1 accepted""")
	parser.add_argument('--read_length', type=int, metavar='INT',
		help="""Trim reads to READ_LENGTH and discard reads with length < READ_LENGTH\nBy default, reads are not trimmed or filtered""")

	args = vars(parser.parse_args())
	return args

def print_species_arguments(args):
	print ("===========Parameters===========")
	print ("Script: run_phylo_cnv.py species")
	print ("Output directory: %s" % args['outdir'])
	print ("Input reads (1st mate): %s" % args['m1'])
	print ("Input reads (2nd mate): %s" % args['m2'])
	print ("Database type: %s" % args['db_type'])
	print ("Remove temporary files: %s" % args['remove_temp'])
	print ("Word size for database search: %s" % args['word_size'])
	if args['mapid']: print ("Minimum mapping identity: %s" % args['mapid'])
	print ("Minimum alignment coverage: %s" % args['aln_cov'])
	print ("Number of reads to use from input: %s" % (args['max_reads'] if args['max_reads'] else 'use all'))
	if args['read_length']: print ("Trim reads to %s-bp and discard reads with length < %s-bp" % (args['read_length'], args['read_length']))
	print ("Number of threads for database search: %s" % args['threads'])

def check_species(args):
	if not os.path.isdir('%s/species' % args['outdir']):
		os.makedirs('%s/species' % args['outdir'])
	if args['word_size'] < 12:
		sys.exit("\nInvalid word size: %s. Must be greater than or equal to 12" % args['word_size'])
	if args['mapid'] and (args['mapid'] < 0 or args['mapid'] > 100):
		sys.exit("\nInvalid mapping identity: %s. Must be between 0 and 100" % args['mapid'])
	if args['aln_cov'] < 0 or args['aln_cov'] > 1:
		sys.exit("\nInvalid alignment coverage: %s. Must be between 0 and 1" % args['aln_cov'])
	for arg in ['m1', 'm2']:
		if args[arg] and not os.path.isfile(args[arg]):
			sys.exit("\nInput file does not exist: '%s'" % args[arg])

def gene_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: run_phylo_cnv.py genes outdir [options]

Description:
This script will map metagenomic reads to bacterial pangenomes and quantify these genes in your data
You can either target one or more specific species, or provide this script a species abundance file
The pipeline can be broken down into three main steps:
  1) build a database of pangenomes for abundance bacterial species
  2) map metagenomic reads to the database
  3) use mapped reads to quantify pangenome genes
""",
		epilog="""Examples:
1) run entire pipeline using defaults:
run_phylo_cnv.py genes /path/to/outdir -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz
			
2) run entire pipeline for a specific species:
run_phylo_cnv.py genes /path/to/outdir --sp_id 57955 -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz

3) just align reads, use faster alignment, only use the first 10M reads, use 4 CPUs:
run_phylo_cnv.py genes /path/to/outdir --align -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz -s very-fast -n 10000000 -t 4

4) just quantify genes, keep reads with >=95% alignment identity and reads with an average quality-score >=30:
run_phylo_cnv.py snps /path/to/outdir --call_genes --mapid 95 --readq 20
	""")

	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str, help='Path to directory to store results. Name should correspond to sample identifier. ')

	parser.add_argument('--tax_mask', action='store_true', default=False, help=argparse.SUPPRESS)
	parser.add_argument('--remove_temp', default=False, action='store_true',
		help='remove intermediate files generated by PhyloCNV. \nintermediate files can be useful to quickly rerun parts of pipeline.')
	
	pipe = parser.add_argument_group('Pipeline options (choose one or more; default=all)')
	pipe.add_argument('--build_db', action='store_true', dest='build_db',
		default=False, help='Build bowtie2 database of pangenomes')
	pipe.add_argument('--align', action='store_true', dest='align',
		default=False, help='Align reads to pangenome database')
	pipe.add_argument('--call_genes', action='store_true', dest='cov',
		default=False, help='Compute coverage of genes in pangenome database')

	db = parser.add_argument_group('Database options (if using --build_db)')
	db.add_argument('-D', type=str, dest='db', help=argparse.SUPPRESS)
	db.add_argument('--sp_cov', type=float, dest='gc_cov', help='Include species with >X coverage (3.0)')
	db.add_argument('--sp_topn', type=int, dest='gc_topn', help='Include top N most abundant species')
	db.add_argument('--sp_id', type=str, dest='gc_id', help='One or more species identifiers to include in database. Separate ids with a comma')
			
	align = parser.add_argument_group('Read alignment options (if using --align)')
	align.add_argument('-1', type=str, dest='m1',
		help='FASTA/FASTQ file containing 1st mate if paired or unpaired reads')
	align.add_argument('-2', type=str, dest='m2',
		help='FASTA/FASTQ file containing 2nd mate if paired')
	align.add_argument('-s', type=str, dest='speed', default='very-sensitive',
		choices=['very-fast', 'fast', 'sensitive', 'very-sensitive'],
		help='Alignment speed/sensitivity (very-sensitive)')
	align.add_argument('-n', type=int, dest='max_reads',
		help='# reads to use from input file(s) (use all)')
	align.add_argument('-t', dest='threads', default=1,
		help='Number of threads to use')
	
	map = parser.add_argument_group('Quantify genes options (if using --call_genes)')
	map.add_argument('--readq', type=int,
		default=20, help='Discard reads with mean quality < READQ (20)')
	map.add_argument('--mapid', type=float,
		default=94.0, help='Discard reads with alignment identity < MAPID (94.0)')
	map.add_argument('--mapq', type=int,
		default=20, help='Discard reads with mapping quality < MAPQ (10)')
	map.add_argument('--aln_cov', type=float,
		default=0.75, help='Discard reads with alignment coverage < ALN_COV (0.75)')
	map.add_argument('--trim', metavar='INT', type=int, default=0,
		help='Trim N base-pairs from read-tails (0)')
		
	args = vars(parser.parse_args())
	if args['gc_id']: args['gc_id'] = args['gc_id'].split(',')
	
	return args


def print_gene_arguments(args):

	print ("===========Parameters===========")

	print ("Script: run_phylo_cnv.py genes")
	print ("Output directory: %s" % args['outdir'])
	print ("Remove temporary files: %s" % args['remove_temp'])
		
	print ("Pipeline options:")
	if args['build_db']:
		print ("  -build bowtie2 database of pangenomes")
	if args['align']:
		print ("  -align reads to bowtie2 pangenome database")
	if args['cov']:
		print ("  -quantify coverage of pangenomes genes")

	if args['build_db']:
		print ("Database options:")
		if args['gc_topn']:
			print ("  -include top %s most abundant species" % args['gc_topn'])
		if args['gc_cov']:
			print ("  -include all species with >=%sX genome coverage" % args['gc_cov'])
		if args['gc_id']:
			print ("  -include specified species id(s): %s" % args['gc_id'])

	if args['align']:
		print ("Read alignment options:")
		print ("  -input reads (1st mate): %s" % args['m1'])
		print ("  -input reads (2nd mate): %s" % args['m2'])
		print ("  -alignment speed/sensitivity: %s" % args['speed'])
		print ("  -number of reads to use from input: %s" % (args['max_reads'] if args['max_reads'] else 'use all'))
		print ("  -number of threads for database search: %s" % args['threads'])

	if args['cov']:
		print ("Gene coverage options:")
		print ("  -minimum alignment percent identity: %s" % args['mapid'])
		print ("  -minimum alignment coverage of reads: %s" % args['aln_cov'])
		print ("  -minimum read quality score: %s" % args['readq'])
		print ("  -minimum mapping quality score: %s" % args['mapq'])
		print ("  -trim %s base-pairs from read-tails" % args['trim'])


def snp_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Usage: run_phylo_cnv.py snps outdir [options]

Description:
This script will map metagenomic reads to bacterial reference genomes and call SNPs in these genomes
You can either target one or more specific species, or provide this script a species abundance file
The pipeline can be broken down into three main steps:
  1) build a database of genome sequences for abundance bacterial species
  2) map metagenomic reads to the database
  3) use mapped reads to call SNPs
""",
		epilog="""Examples:
1) run entire pipeline using defaults:
run_phylo_cnv.py snps /path/to/outdir -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz
			
2) run entire pipeline for a specific species:
run_phylo_cnv.py snps /path/to/outdir --sp_id 57955 -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz

3) just align reads, use faster alignment, only use the first 10M reads, use 4 CPUs:
run_phylo_cnv.py snps /path/to/outdir --align -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz -s very-fast -n 10000000 -t 4

4) just call SNPs, keep reads with >=95% alignment identity and keep bases with quality-scores >=35:
run_phylo_cnv.py snps /path/to/outdir --call_snps --mapid 95 --baseq 35
	""")
	
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str, help='Path to directory to store results. Name should correspond to sample identifier.')
	
	parser.add_argument('--tax_mask', action='store_true', default=False, help=argparse.SUPPRESS)
	parser.add_argument('--remove_temp', default=False, action='store_true',
		help='remove intermediate files generated by PhyloCNV. \nintermediate files can be useful to quickly rerun parts of pipeline.')
	
	pipe = parser.add_argument_group('Pipeline options (choose one or more; default=all)')
	pipe.add_argument('--build_db', action='store_true', dest='build_db',
		default=False, help='Build bowtie2 database of pangenomes')
	pipe.add_argument('--align', action='store_true', dest='align',
		default=False, help='Align reads to pangenome database')
	pipe.add_argument('--call_snps', action='store_true', dest='call',
		default=False, help='Run samtools mpileup and call SNPs')
				
	db = parser.add_argument_group('Database options (if using --build_db)')
	db.add_argument('-D', type=str, dest='db', help=argparse.SUPPRESS)
	db.add_argument('--sp_cov', type=float, dest='gc_cov', help='Include species with >X coverage (3.0)')
	db.add_argument('--sp_topn', type=int, dest='gc_topn', help='Include top N most abundant species')
	db.add_argument('--sp_id', type=str, dest='gc_id', help='One or more species identifiers to include in database. Separate ids with a comma')
			
	align = parser.add_argument_group('Read alignment options (if using --align)')
	align.add_argument('-1', type=str, dest='m1', help='FASTA/FASTQ file containing 1st mate if paired or unpaired reads')
	align.add_argument('-2', type=str, dest='m2', help='FASTA/FASTQ file containing 2nd mate if paired')
	align.add_argument('-s', type=str, dest='speed', default='very-sensitive',
		choices=['very-fast', 'fast', 'sensitive', 'very-sensitive'],
		help='Bowtie2 alignment speed/sensitivity (very-sensitive)')
	align.add_argument('-n', type=int, dest='max_reads', help='# reads to use from input file(s) (use all)')
	align.add_argument('-t', dest='threads', default=1, help='Number of threads to use')
	
	snps = parser.add_argument_group('SNP calling options (if using --call_snps)')
	snps.add_argument('--mapid', type=float,
		default=94.0, help='Discard reads with alignment identity < MAPID (94.0)')
	snps.add_argument('--mapq', type=int,
		default=20, help='Discard reads with mapping quality < MAPQ (20)')
	snps.add_argument('--baseq', type=int,
		default=30, help='Discard bases with quality < BASEQ (30)')
	snps.add_argument('--readq', type=int,
		default=20, help='Discard reads with mean quality < READQ (20)')
	snps.add_argument('--trim', metavar='INT', type=int, default=0,
		help='Trim N base-pairs from read-tails (0)')
	snps.add_argument('--baq', default=False, action='store_true',
		help='Enable BAQ (per-base alignment quality)')
	snps.add_argument('--redo_baq', default=False, action='store_true',
		help='Recalculate BAQ on the fly')
	snps.add_argument('--adjust_mq', default=False, action='store_true',
		help='Adjust MAPQ')
			
	args = vars(parser.parse_args())
	if args['gc_id']: args['gc_id'] = args['gc_id'].split(',')
	
	return args

def print_snp_arguments(args):

	print ("===========Parameters===========")

	print ("Script: run_phylo_cnv.py snps")
	print ("Output directory: %s" % args['outdir'])
	print ("Remove temporary files: %s" % args['remove_temp'])
		
	print ("Pipeline options:")
	if args['build_db']:
		print ("  -build bowtie2 database of genomes")
	if args['align']:
		print ("  -align reads to bowtie2 genome database")
	if args['call']:
		print ("  -use samtools to generate pileups and call SNPs")

	if args['build_db']:
		print ("Database options:")
		if args['gc_topn']:
			print ("  -include top %s most abundant species" % args['gc_topn'])
		if args['gc_cov']:
			print ("  -include all species with >=%sX genome coverage" % args['gc_cov'])
		if args['gc_id']:
			print ("  -include specified species id(s): %s" % args['gc_id'])

	if args['align']:
		print ("Read alignment options:")
		print ("  -input reads (1st mate): %s" % args['m1'])
		print ("  -input reads (2nd mate): %s" % args['m2'])
		print ("  -alignment speed/sensitivity: %s" % args['speed'])
		print ("  -number of reads to use from input: %s" % (args['max_reads'] if args['max_reads'] else 'use all'))
		print ("  -number of threads for database search: %s" % args['threads'])

	if args['call']:
		print ("SNP calling options:")
		print ("  -minimum alignment percent identity: %s" % args['mapid'])
		print ("  -minimum mapping quality score: %s" % args['mapq'])
		print ("  -minimum base quality score: %s" % args['baseq'])
		print ("  -minimum read quality score: %s" % args['readq'])
		print ("  -trim %s base-pairs from read-tails" % args['trim'])
		if args['baq']: print ("  -enable BAQ (per-base alignment quality)")
		if args['redo_baq']: print ("  -recalculate BAQ on the fly")
		if args['adjust_mq']: print ("  -adjust MAPQ")

def check_genes(args):
	""" Check validity of command line arguments """
	# create output directory
	if not os.path.isdir('%s/genes' % args['outdir']):
		os.makedirs('%s/genes' % args['outdir'])
	# turn on pipeline options
	if not any([args['build_db'], args['align'], args['cov']]):
		args['build_db'] = True
		args['align'] = True
		args['cov'] = True
	# set default species selection
	if not any([args['gc_id'], args['gc_topn'], args['gc_cov']]):
		args['gc_cov'] = 3.0
	# species selection options, but no no profile file
	profile='%s/species/species_profile.txt' % args['outdir']
	if not os.path.isfile(profile):
		if (args['gc_topn'] or args['gc_cov']) and args['build_db']:
			sys.exit("\nCould not find species abundance profile: %s\nTo specify species with --sp_topn or --sp_cov you must have run: run_phylo_cnv.py species" % profile)
	# no database but --align specified
	if (args['align']
		and not args['build_db']
		and not os.path.isfile('%s/genes/db/pangenomes.fa' % args['outdir'])):
		error = "\nYou've specified --align, but no database has been built"
		error += "\nTry running with --build_db"
		sys.exit(error)
	# no bamfile but --cov specified
	if (args['cov']
		and not args['align']
		and not os.path.isfile('%s/genes/pangenome.bam' % args['outdir'])):
		error = "\nYou've specified --call_genes, but no alignments were found"
		error += "\nTry running with --align"
		sys.exit(error)
	# no reads
	if args['align'] and not args['m1']:
		sys.exit("\nTo align reads, you must specify path to input FASTA/FASTQ")
	# check input file paths
	for arg in ['m1', 'm2']:
		if args[arg] and not os.path.isfile(args[arg]):
			sys.exit("\nInput file does not exist: '%s'" % args[arg])
	# input options
	if args['m2'] and not args['m1']:
		sys.exit("\nMust specify -1 and -2 if aligning paired end reads")
	# sanity check input values
	if args['mapid'] < 1 or args['mapid'] > 100:
		sys.exit("\nMAPID must be between 1 and 100")
	if args['aln_cov'] < 0 or args['aln_cov'] > 1:
		sys.exit("\nALN_COV must be between 0 and 1")

def check_snps(args):
	""" Check validity of command line arguments """
	# create output directory
	if not os.path.isdir('%s/snps' % args['outdir']):
		os.makedirs('%s/snps' % args['outdir'])
	# pipeline options
	if not any([args['build_db'], args['align'], args['call']]):
		args['build_db'] = True
		args['align'] = True
		args['call'] = True
	# set default species selection
	if not any([args['gc_id'], args['gc_topn'], args['gc_cov']]):
		args['gc_cov'] = 3.0
	# species selection options, but no no profile file
	profile='%s/species/species_profile.txt' % args['outdir']
	if not os.path.isfile(profile):
		if (args['gc_topn'] or args['gc_cov']) and args['build_db']:
			sys.exit("\nCould not find species abundance profile: %s\nTo specify species with --sp_topn or --sp_cov you must have run: run_phylo_cnv.py species" % profile)
	# no database but --align specified
	if (args['align']
		and not args['build_db']
		and not os.path.isfile('%s/snps/db/genomes.fa' % args['outdir'])):
		error = "\nYou've specified --align, but no database has been built"
		error += "\nTry running with --build_db"
		sys.exit(error)
	# no bamfile but --call specified
	if (args['call']
		and not args['align']
		and not os.path.isfile('%s/snps/genomes.bam' % args['outdir'])
		):
		error = "\nYou've specified --call_snps, but no alignments were found"
		error += "\nTry running with --align"
		sys.exit(error)
	# no genomes but --call specified
	if (args['call']
		and not args['build_db']
		and not os.path.isfile('%s/snps/db/genomes.fa' % args['outdir'])
		):
		error = "\nYou've specified --call_snps, but the no genome database was found"
		error += "\nTry running with --build_db"
		sys.exit(error)
	# no reads
	if args['align'] and not args['m1']:
		sys.exit("\nTo align reads, you must specify path to input FASTA/FASTQ")
	# check input file paths
	for arg in ['m1', 'm2']:
		if args[arg] and not os.path.isfile(args[arg]):
			sys.exit("\nInput file does not exist: '%s'" % args[arg])
	# input options
	if args['m2'] and not args['m1']:
		sys.exit("\nMust specify -1 and -2 if aligning paired end reads")
	# sanity check input values
	if args['mapid'] < 1 or args['mapid'] > 100:
		sys.exit("\nMAPQ must be between 1 and 100")
	if args['mapq'] < 0 or args['mapq'] > 100:
		sys.exit("\nMAPQ must be between 0 and 100")
	if args['baseq'] < 0 or args['baseq'] > 100:
		sys.exit("\nBASEQ must be between 0 and 100")

if __name__ == '__main__':

	program = get_program()
	args = get_arguments(program)
	check_arguments(program, args)
	utility.print_copyright()
	print_arguments(program, args)
	run_program(program, args)
	



