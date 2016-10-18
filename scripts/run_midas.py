#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import argparse, sys, os, platform
from midas import utility

def get_program():
	""" Get program specified by user (species, genes, or snps) """
	if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
		print('')
		print('Usage: run_midas.py <command> [options]')
		print('')
		print('Note: use run_midas.py <command> -h to view usage for a specific command')
		print('')
		print('Commands:')
		print('\tspecies\t estimate the abundance of 5,952 bacterial species')
		print('\tgenes\t identify gene copy number variants in abundant species')
		print('\tsnps\t identify single nucleotide variants in abundant species')
		quit()
	elif sys.argv[1] not in ['species', 'genes', 'snps']:
		sys.exit("Error: Unrecognized command: '%s'" % sys.argv[1])
		quit()
	else:
		return sys.argv[1]

def open_log(program, args):
	""" Open log file """
	logpath = '%s/%s/%s' % (args['outdir'], program, 'log.txt')
	args['log'] = open(logpath, 'w')

def get_arguments(program):
	""" Get arguments for specified program """
	if program == 'species':
		args = species_arguments()
	elif program == 'genes':
		args = gene_arguments()
	elif program == 'snps':
		args = snp_arguments()
	else:
		sys.exit("Error: Unrecognized program: '%s'" % program)
	utility.add_executables(args)
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
		sys.exit("Error: Unrecognized program: '%s'" % program)
	if platform.system() not in ['Linux', 'Darwin']:
		sys.exit("Error: Operating system '%s' not supported" % system())

def print_arguments(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		print_species_arguments(args)
	elif program == 'genes':
		print_gene_arguments(args)
	elif program == 'snps':
		print_snp_arguments(args)
	else:
		sys.exit("Error: Unrecognized program: '%s'" % program)
	
def run_program(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		from midas.run import species
		species.run_pipeline(args)
	elif program == 'genes':
		from midas.run import genes
		genes.run_pipeline(args)
	elif program == 'snps':
		from midas.run import snps
		snps.run_pipeline(args)
	else:
		sys.exit("Error: Unrecognized program: '%s'" % program)

def species_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description:
This script will map metagenomic reads to a database of phylogenetic marker genes using HS-BLASTN
Mapped reads are used estimate the read depth and relative abundance of bacterial species
Reads are mapped according to gene-specific, species-level mapping thresholds (94.5-98% DNA identity)
Reads that map equally well to 2 or more species are probabalistically assigned

Usage: run_midas.py species outdir [options]
""",
		epilog="""Examples:
1) run with defaults using a paired-end metagenome:
run_midas.py species /path/to/outdir -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz

2) run using a single-end metagenome with 4 CPUs and only 4M reads:
run_midas.py species /path/to/outdir -1 /path/to/reads_1.fq.gz -t 4 -n 4000000

3) run with exactly 80 base-pair reads:
run_midas.py species /path/to/outdir -1 /path/to/reads_1.fq.gz --read_length 80

4) quantify species abundance using a 16S database:
run_midas.py species /path/to/outdir -1 /path/to/reads_1.fq.gz --db_type ssuRNA

Output files:
1) species_profile.txt: tab-delimited output file containing abundances of 5,952 species
2) log.txt: log file containing parameters used
3) temp/: directory of intermediate files. use `--remove_temp` to remove these files

Output formats:
species_profile.txt
1) species_id: species (i.e. genome-cluster) identifier
2) species_name: unique species name
3) count_reads: number of reads mapped to marker genes
4) coverage: estimated genome-coverage of species in metagenome
5) relative_abundance: estimated relative abundance of species in metagenome
""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str, help='Path to directory to store results. Name should correspond to sample identifier.')
	parser.add_argument('-1', type=str, dest='m1', required=True,
		help="""FASTA/FASTQ file containing 1st mate if using paired-end reads.
Otherwise FASTA/FASTQ containing unpaired reads.
Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)""")
	parser.add_argument('-2', type=str, dest='m2',
		help="""FASTA/FASTQ file containing 2nd mate if using paired-end reads.
Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)""")
	parser.add_argument('-n', type=int, dest='max_reads',
		help="""Number of reads to use from input file(s) (use all)""")
	parser.add_argument('-t', dest='threads', default=1,
		help="""Number of threads to use for database search (1)""")
	parser.add_argument('-d', type=str, dest='db', default=os.environ['MIDAS_DB'] if 'MIDAS_DB' in os.environ else None,
		help="""Path to reference database""")
	parser.add_argument('--db_type', choices=['phyeco', 'ssuRNA'], metavar='', default='phyeco',
		help=argparse.SUPPRESS)
		#help="""Reference database. Choices:\n'phyeco': universal-single-copy protein family database (default)\n'ssuRNA': 16S ribosomal rna database""")
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
	lines = []
	lines.append("===========Parameters===========")
	lines.append("Script: run_midas.py species")
	lines.append("Output directory: %s" % args['outdir'])
	lines.append("Input reads (1st mate): %s" % args['m1'])
	lines.append("Input reads (2nd mate): %s" % args['m2'])
	lines.append("Database type: %s" % args['db_type'])
	lines.append("Remove temporary files: %s" % args['remove_temp'])
	lines.append("Word size for database search: %s" % args['word_size'])
	if args['mapid']: lines.append("Minimum mapping identity: %s" % args['mapid'])
	lines.append("Minimum alignment coverage: %s" % args['aln_cov'])
	lines.append("Number of reads to use from input: %s" % (args['max_reads'] if args['max_reads'] else 'use all'))
	if args['read_length']: lines.append("Trim reads to %s-bp and discard reads with length < %s-bp" % (args['read_length'], args['read_length']))
	lines.append("Number of threads for database search: %s" % args['threads'])
	args['log'].write('\n'.join(lines)+'\n')
	sys.stdout.write('\n'.join(lines)+'\n')

def check_species(args):
	# check file type
	if args['m1']: args['file_type'] = utility.auto_detect_file_type(args['m1'])
	# check database
	check_database(args)
	# create output directories
	if not os.path.isdir('%s/species' % args['outdir']):
		os.makedirs('%s/species' % args['outdir'])
	# check word size
	if args['word_size'] < 12:
		sys.exit("\nError: Invalid word size: %s. Must be greater than or equal to 12" % args['word_size'])
	# check mapping identity
	if args['mapid'] and (args['mapid'] < 0 or args['mapid'] > 100):
		sys.exit("\nError: Invalid mapping identity: %s. Must be between 0 and 100" % args['mapid'])
	# check alignment coverage
	if args['aln_cov'] < 0 or args['aln_cov'] > 1:
		sys.exit("\nError: Invalid alignment coverage: %s. Must be between 0 and 1" % args['aln_cov'])
	# check that m1 (and m2) exist
	for arg in ['m1', 'm2']:
		if args[arg] and not os.path.isfile(args[arg]):
			sys.exit("\nError: Input file does not exist: '%s'" % args[arg])
	# check that extention matches compression
	if args['m1']: utility.check_compression(args['m1'])
	if args['m2']: utility.check_compression(args['m2'])

def create_directories(program, args):
	dirs = [args['outdir']]
	dirs.append('%s/%s' % (args['outdir'], program))
	if program != 'species':
		dirs.append('%s/%s/%s' % (args['outdir'], program, 'output'))
	dirs.append('%s/%s/%s' % (args['outdir'], program, 'temp'))
	for dir in dirs:
		if not os.path.isdir(dir): os.mkdir(dir)

def gene_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description:
This script will map metagenomic reads to bacterial pangenomes and quantify these genes in your data
The pipeline can be broken down into three main steps:
  1) build a database of pangenomes for abundance bacterial species
  2) map metagenomic reads to the database
  3) use mapped reads to quantify pangenome genes
  
Usage: run_midas.py genes outdir [options]
""",
		epilog="""Examples:
1) run entire pipeline using defaults:
run_midas.py genes /path/to/outdir -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz
			
2) run entire pipeline for a specific species:
run_midas.py genes /path/to/outdir --species_id 57955 -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz

3) just align reads, use faster alignment, only use the first 10M reads, use 4 CPUs:
run_midas.py genes /path/to/outdir --align -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz -s very-fast -n 10000000 -t 4

4) just quantify genes, keep reads with >=95% alignment identity and reads with an average quality-score >=30:
run_midas.py snps /path/to/outdir --call_genes --mapid 95 --readq 20

Output files:
1) output/: directory of per-species output files. are tab-delimited and gzip-compressed. named with the convention {SPECIES_ID}.genes.gz
2) summary.txt: tab-delimited file summarizing alignments
3) log.txt: log file containing parameters used
4) temp/: directory of intermediate files. use `--remove_temp` to remove these files

Output formats:
output/
1) gene_id: identifer of gene family (99% identity clustering cutoff); 'peg' and 'rna' indicate coding & non-coding genes respectively
2) coverage: average read-depth of gene_id (# aligned bp / gene length in bp)
3) copy_number: estimated copy-number of gene_id in population  (coverage of gene_id / median coverage of 15 universal single copy genes)

summary.txt
1) species_id: species identifer
2) pangenome_size: total number of gene families (99% identity clustering cutoff) in reference pangenome
3) covered_genes: number of pangenome gene families with non-zero coverage
4) fraction_covered: fraction of pangenome gene families with non-zero coverage
5) mean_coverage: mean read-depth across gene families with non-zero coverage
6) marker_coverage: median read-depth across 15 universal single copy genes
""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str, help='Path to directory to store results. Name should correspond to sample identifier. ')
	parser.add_argument('--remove_temp', default=False, action='store_true',
		help='remove intermediate files generated by MIDAS. \nintermediate files can be useful to quickly rerun parts of pipeline.')
	pipe = parser.add_argument_group('Pipeline options (choose one or more; default=all)')
	pipe.add_argument('--build_db', action='store_true', dest='build_db',
		default=False, help='Build bowtie2 database of pangenomes')
	pipe.add_argument('--align', action='store_true', dest='align',
		default=False, help='Align reads to pangenome database')
	pipe.add_argument('--call_genes', action='store_true', dest='cov',
		default=False, help='Compute coverage of genes in pangenome database')
	db = parser.add_argument_group('Database options (if using --build_db)')
	db.add_argument('-d', type=str, dest='db', default=os.environ['MIDAS_DB'] if 'MIDAS_DB' in os.environ else None,
		help="""Path to reference database""")
	db.add_argument('--species_cov', type=float, dest='species_cov', metavar='FLOAT', help='Include species with >X coverage (3.0)')
	db.add_argument('--species_topn', type=int, dest='species_topn', metavar='INT', help='Include top N most abundant species')
	db.add_argument('--species_id', type=str, dest='species_id', metavar='CHAR', help='One or more species identifiers to include in database. Separate ids with a comma')
	align = parser.add_argument_group('Read alignment options (if using --align)')
	align.add_argument('-1', type=str, dest='m1', required=True,
		help="""FASTA/FASTQ file containing 1st mate if using paired-end reads.
Otherwise FASTA/FASTQ containing unpaired reads.
Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)""")
	align.add_argument('-2', type=str, dest='m2',
		help="""FASTA/FASTQ file containing 2nd mate if using paired-end reads.
Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)""")
	align.add_argument('-s', type=str, dest='speed', default='very-sensitive',
		choices=['very-fast', 'fast', 'sensitive', 'very-sensitive'],
		help='Alignment speed/sensitivity (very-sensitive)')
	align.add_argument('-n', type=int, dest='max_reads',
		help='# reads to use from input file(s) (use all)')
	align.add_argument('-t', dest='threads', default=1,
		help='Number of threads to use')
	map = parser.add_argument_group('Quantify genes options (if using --call_genes)')
	map.add_argument('--readq', type=int, metavar='INT',
		default=20, help='Discard reads with mean quality < READQ (20)')
	map.add_argument('--mapid', type=float, metavar='FLOAT',
		default=94.0, help='Discard reads with alignment identity < MAPID (94.0)')
	map.add_argument('--mapq', type=int, metavar='INT',
		default=20, help='Discard reads with mapping quality < MAPQ (10)')
	map.add_argument('--aln_cov', type=float, metavar='FLOAT',
		default=0.75, help='Discard reads with alignment coverage < ALN_COV (0.75)')
	map.add_argument('--trim', type=int, default=0, metavar='INT',
		help='Trim N base-pairs from read-tails (0)')
	args = vars(parser.parse_args())
	if args['species_id']: args['species_id'] = args['species_id'].split(',')
	return args

def print_gene_arguments(args):
	lines = []
	lines.append("===========Parameters===========")
	lines.append("Script: run_midas.py genes")
	lines.append("Output directory: %s" % args['outdir'])
	lines.append("Remove temporary files: %s" % args['remove_temp'])
	lines.append("Pipeline options:")
	if args['build_db']:
		lines.append("  -build bowtie2 database of pangenomes")
	if args['align']:
		lines.append("  -align reads to bowtie2 pangenome database")
	if args['cov']:
		lines.append("  -quantify coverage of pangenomes genes")
	if args['build_db']:
		lines.append("Database options:")
		if args['species_topn']:
			lines.append("  -include top %s most abundant species" % args['species_topn'])
		if args['species_cov']:
			lines.append("  -include all species with >=%sX genome coverage" % args['species_cov'])
		if args['species_id']:
			lines.append("  -include specified species id(s): %s" % args['species_id'])
	if args['align']:
		lines.append("Read alignment options:")
		lines.append("  -input reads (1st mate): %s" % args['m1'])
		lines.append("  -input reads (2nd mate): %s" % args['m2'])
		lines.append("  -alignment speed/sensitivity: %s" % args['speed'])
		lines.append("  -number of reads to use from input: %s" % (args['max_reads'] if args['max_reads'] else 'use all'))
		lines.append("  -number of threads for database search: %s" % args['threads'])
	if args['cov']:
		lines.append("Gene coverage options:")
		lines.append("  -minimum alignment percent identity: %s" % args['mapid'])
		lines.append("  -minimum alignment coverage of reads: %s" % args['aln_cov'])
		lines.append("  -minimum read quality score: %s" % args['readq'])
		lines.append("  -minimum mapping quality score: %s" % args['mapq'])
		lines.append("  -trim %s base-pairs from read-tails" % args['trim'])
	args['log'].write('\n'.join(lines)+'\n')
	sys.stdout.write('\n'.join(lines)+'\n')

def snp_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description:
This script will map metagenomic reads to bacterial reference genomes and call SNPs in these genomes
The pipeline can be broken down into three main steps:
  1) build a database of genome sequences for abundance bacterial species
  2) map metagenomic reads to the database
  3) use mapped reads to call SNPs
  
Usage: run_midas.py snps outdir [options]
""",
		epilog="""Examples:
1) run entire pipeline using defaults:
run_midas.py snps /path/to/outdir -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz
			
2) run entire pipeline for a specific species:
run_midas.py snps /path/to/outdir --species_id 57955 -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz

3) just align reads, use faster alignment, only use the first 10M reads, use 4 CPUs:
run_midas.py snps /path/to/outdir --align -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz -s very-fast -n 10000000 -t 4

4) just call SNPs, keep reads with >=95% alignment identity and keep bases with quality-scores >=35:
run_midas.py snps /path/to/outdir --call_snps --mapid 95 --baseq 35
	
Output files:
1) output/: directory of per-species output files. are tab-delimited and gzip-compressed. named with the convention {SPECIES_ID}.snps.gz
2) summary.txt: tab-delimited file summarizing alignments
3) log.txt: log file containing parameters used
4) temp/: directory ofintermediate files. use `--remove_temp` to remove these files

Output formats:
output/
1) ref_id: scaffold id
2) ref_pos: position on scaffold
3) ref_allele: reference allele
4) alt_allele: alternate allele
5) cons_allele: consensus allele
6) count_alleles: number of distinct alleles observed
7) count_ref: count reference alleles observed
8) count_alt: count alternate alleles observed
9) depth: count total reads at ref_pos
10) ref_freq: frequency (0.0 to 1.0) of reference allele

summary.txt:
1) species_id: species identifer
2) genome_length: length of reference genome used for read-mapping
3) covered_bases: number of genomic positions covered
4) fraction_covered: fraction of genomic positions with non-zero coverage
5) mean_coverage: mean read-depth at covered genomic positions
""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str, help='Path to directory to store results. Name should correspond to sample identifier.')
	parser.add_argument('--remove_temp', default=False, action='store_true',
		help='remove intermediate files generated by MIDAS. \nintermediate files can be useful to quickly rerun parts of pipeline.')
	pipe = parser.add_argument_group('Pipeline options (choose one or more; default=all)')
	pipe.add_argument('--build_db', action='store_true', dest='build_db',
		default=False, help='Build bowtie2 database of pangenomes')
	pipe.add_argument('--align', action='store_true', dest='align',
		default=False, help='Align reads to pangenome database')
	pipe.add_argument('--call_snps', action='store_true', dest='call',
		default=False, help='Run samtools mpileup and call SNPs')
	db = parser.add_argument_group('Database options (if using --build_db)')
	db.add_argument('-d', type=str, dest='db', default=os.environ['MIDAS_DB'] if 'MIDAS_DB' in os.environ else None,
		help="""Path to reference database""")
	db.add_argument('--species_cov', type=float, dest='species_cov', metavar='FLOAT', help='Include species with >X coverage (3.0)')
	db.add_argument('--species_topn', type=int, dest='species_topn', metavar='INT', help='Include top N most abundant species')
	db.add_argument('--species_id', type=str, dest='species_id', metavar='CHAR', help='One or more species identifiers to include in database. Separate ids with a comma')
	align = parser.add_argument_group('Read alignment options (if using --align)')
	align.add_argument('-1', type=str, dest='m1', required=True,
		help="""FASTA/FASTQ file containing 1st mate if using paired-end reads.
Otherwise FASTA/FASTQ containing unpaired reads.
Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)""")
	align.add_argument('-2', type=str, dest='m2',
		help="""FASTA/FASTQ file containing 2nd mate if using paired-end reads.
Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)""")
	align.add_argument('-s', type=str, dest='speed', default='very-sensitive',
		choices=['very-fast', 'fast', 'sensitive', 'very-sensitive'],
		help='Bowtie2 alignment speed/sensitivity (very-sensitive)')
	align.add_argument('-n', type=int, dest='max_reads', help='# reads to use from input file(s) (use all)')
	align.add_argument('-t', dest='threads', default=1, help='Number of threads to use')
	snps = parser.add_argument_group('SNP calling options (if using --call_snps)')
	snps.add_argument('--mapid', type=float, metavar='FLOAT',
		default=94.0, help='Discard reads with alignment identity < MAPID (94.0)')
	snps.add_argument('--mapq', type=int, metavar='INT',
		default=20, help='Discard reads with mapping quality < MAPQ (20)')
	snps.add_argument('--baseq', type=int, metavar='INT',
		default=30, help='Discard bases with quality < BASEQ (30)')
	snps.add_argument('--readq', type=int, metavar='INT',
		default=20, help='Discard reads with mean quality < READQ (20)')
	snps.add_argument('--trim', metavar='INT', type=int, default=0,
		help='Trim N base-pairs from read-tails (0)')
	snps.add_argument('--discard', default=False, action='store_true',
		help='Discard discordant read-pairs')
	snps.add_argument('--baq', default=False, action='store_true',
		help='Enable BAQ (per-base alignment quality)')
	snps.add_argument('--adjust_mq', default=False, action='store_true',
		help='Adjust MAPQ')
	args = vars(parser.parse_args())
	if args['species_id']: args['species_id'] = args['species_id'].split(',')
	return args

def print_snp_arguments(args):
	lines = []
	lines.append("===========Parameters===========")
	lines.append("Script: run_midas.py snps")
	lines.append("Output directory: %s" % args['outdir'])
	lines.append("Remove temporary files: %s" % args['remove_temp'])
	lines.append("Pipeline options:")
	if args['build_db']:
		lines.append("  -build bowtie2 database of genomes")
	if args['align']:
		lines.append("  -align reads to bowtie2 genome database")
	if args['call']:
		lines.append("  -use samtools to generate pileups and call SNPs")
	if args['build_db']:
		lines.append("Database options:")
		if args['species_topn']:
			lines.append("  -include top %s most abundant species" % args['species_topn'])
		if args['species_cov']:
			lines.append("  -include all species with >=%sX genome coverage" % args['species_cov'])
		if args['species_id']:
			lines.append("  -include specified species id(s): %s" % args['species_id'])
	if args['align']:
		lines.append("Read alignment options:")
		lines.append("  -input reads (1st mate): %s" % args['m1'])
		lines.append("  -input reads (2nd mate): %s" % args['m2'])
		lines.append("  -alignment speed/sensitivity: %s" % args['speed'])
		lines.append("  -number of reads to use from input: %s" % (args['max_reads'] if args['max_reads'] else 'use all'))
		lines.append("  -number of threads for database search: %s" % args['threads'])
	if args['call']:
		lines.append("SNP calling options:")
		lines.append("  -minimum alignment percent identity: %s" % args['mapid'])
		lines.append("  -minimum mapping quality score: %s" % args['mapq'])
		lines.append("  -minimum base quality score: %s" % args['baseq'])
		lines.append("  -minimum read quality score: %s" % args['readq'])
		lines.append("  -trim %s base-pairs from read-tails" % args['trim'])
		if args['discard']: lines.append("  -discard discordant read-pairs")
		if args['baq']: lines.append("  -enable BAQ (per-base alignment quality)")
		if args['adjust_mq']: lines.append("  -adjust MAPQ")
	args['log'].write('\n'.join(lines)+'\n')
	sys.stdout.write('\n'.join(lines)+'\n')

def check_database(args):
	if 'db' is None:
		error = "\nError: No reference database specified\n"
		error = "Use the flag -d to specify a database,\n"
		error = "or set the MIDAS_DB environmental variable: export MIDAS_DB=/path/to/midas/db\n"
		sys.exit(error)
	if not os.path.isdir(args['db']):
		error = "\nError: Specified reference database does not exist: %s\n" % args['db']
		error += "\nCheck that you've entered the path correctly and the database exists"
		error += "\nTo download the default database, run: MIDAS/scripts/download_ref_db.py"
		error += "\nTo build a custom database, run: MIDAS/scripts/build_midas_db.py\n"
		sys.exit(error)
	for dir in ['marker_genes', 'pan_genomes', 'rep_genomes']:
		path = '%s/%s' % (args['db'], dir)
		if not os.path.isdir(path):
			error = "\nError: Could not locate required database directory: %s\n" % path
			sys.exit(error)
	for file in ['species_info.txt']:
		path = '%s/%s' % (args['db'], file)
		if not os.path.exists(path):
			error = "\nError: Could not locate required database file: %s\n" % path
			sys.exit(error)

def check_selected_species(args):
	if args['species_id']:
		for species_id in args['species_id']:
			if args['program'] == 'genes':
				path = '%s/pan_genomes/%s' % (args['db'], species_id)
				error = "\nError: Could not locate species pan genome: %s\n" % path
			if args['program'] == 'snps':
				path = '%s/rep_genomes/%s' % (args['db'], species_id)
				error = "\nError: Could not locate species representative genome: %s\n" % path
			if not os.path.isdir(path):
				 sys.exit(error)

def check_genes(args):
	""" Check validity of command line arguments """
	# check file type
	if args['m1']: args['file_type'] = utility.auto_detect_file_type(args['m1'])
	# check database
	check_database(args)
	# make sure selected species are valid
	check_selected_species(args)
	# create output directory
	if not os.path.isdir('%s/genes' % args['outdir']):
		os.makedirs('%s/genes' % args['outdir'])
	# turn on pipeline options
	if not any([args['build_db'], args['align'], args['cov']]):
		args['build_db'] = True
		args['align'] = True
		args['cov'] = True
	# set default species selection
	if not any([args['species_id'], args['species_topn'], args['species_cov']]):
		args['species_cov'] = 3.0
	# species selection options, but no no profile file
	profile='%s/species/species_profile.txt' % args['outdir']
	if not os.path.isfile(profile):
		if (args['species_topn'] or args['species_cov']) and args['build_db']:
			sys.exit("\nError: Could not find species abundance profile: %s\n\
To specify species with --species_topn or --species_cov you must have run: run_midas.py species\n\
Alternatively, you can manually specify one or more species using --species_id" % profile)
	# no database but --align specified
	if (args['align']
		and not args['build_db']
		and not os.path.isfile('%s/genes/temp/pangenomes.fa' % args['outdir'])):
		error = "\nError: You've specified --align, but no database has been built"
		error += "\nTry running with --build_db"
		sys.exit(error)
	# no bamfile but --cov specified
	if (args['cov']
		and not args['align']
		and not os.path.isfile('%s/genes/temp/pangenome.bam' % args['outdir'])):
		error = "\nError: You've specified --call_genes, but no alignments were found"
		error += "\nTry running with --align"
		sys.exit(error)
	# no reads
	if args['align'] and not args['m1']:
		sys.exit("\nError: To align reads, you must specify path to input FASTA/FASTQ")
	# check input file paths
	for arg in ['m1', 'm2']:
		if args[arg] and not os.path.isfile(args[arg]):
			sys.exit("\nError: Input file does not exist: '%s'" % args[arg])
	# check compression
	if args['m1']: utility.check_compression(args['m1'])
	if args['m2']: utility.check_compression(args['m2'])
	# input options
	if args['m2'] and not args['m1']:
		sys.exit("\nError: Must specify -1 and -2 if aligning paired end reads")
	# sanity check input values
	if args['mapid'] < 1 or args['mapid'] > 100:
		sys.exit("\nError: MAPID must be between 1 and 100")
	if args['aln_cov'] < 0 or args['aln_cov'] > 1:
		sys.exit("\nError: ALN_COV must be between 0 and 1")

def check_snps(args):
	""" Check validity of command line arguments """
	# check file type
	if args['m1']: args['file_type'] = utility.auto_detect_file_type(args['m1'])
	# check database
	check_database(args)
	# make sure selected species are valid
	check_selected_species(args)
	# create output directory
	if not os.path.isdir('%s/snps' % args['outdir']):
		os.makedirs('%s/snps' % args['outdir'])
	# pipeline options
	if not any([args['build_db'], args['align'], args['call']]):
		args['build_db'] = True
		args['align'] = True
		args['call'] = True
	# set default species selection
	if not any([args['species_id'], args['species_topn'], args['species_cov']]):
		args['species_cov'] = 3.0
	# species selection options, but no no profile file
	profile='%s/species/species_profile.txt' % args['outdir']
	if not os.path.isfile(profile):
		if (args['species_topn'] or args['species_cov']) and args['build_db']:
			sys.exit("\nError: Could not find species abundance profile: %s\n\
To specify species with --species_topn or --species_cov you must have run: run_midas.py species\n\
Alternatively, you can manually specify one or more species using --species_id" % profile)
	# no database but --align specified
	if (args['align']
		and not args['build_db']
		and not os.path.isfile('%s/snps/temp/genomes.fa' % args['outdir'])):
		error = "\nError: You've specified --align, but no database has been built"
		error += "\nTry running with --build_db"
		sys.exit(error)
	# no bamfile but --call specified
	if (args['call']
		and not args['align']
		and not os.path.isfile('%s/snps/temp/genomes.bam' % args['outdir'])
		):
		error = "\nError: You've specified --call_snps, but no alignments were found"
		error += "\nTry running with --align"
		sys.exit(error)
	# no genomes but --call specified
	if (args['call']
		and not args['build_db']
		and not os.path.isfile('%s/snps/temp/genomes.fa' % args['outdir'])
		):
		error = "\nError: You've specified --call_snps, but the no genome database was found"
		error += "\nTry running with --build_db"
		sys.exit(error)
	# no reads
	if args['align'] and not args['m1']:
		sys.exit("\nError: To align reads, you must specify path to input FASTA/FASTQ")
	# check input file paths
	for arg in ['m1', 'm2']:
		if args[arg] and not os.path.isfile(args[arg]):
			sys.exit("\nError: Input file does not exist: '%s'" % args[arg])
	# check compression
	if args['m1']: utility.check_compression(args['m1'])
	if args['m2']: utility.check_compression(args['m2'])
	# input options
	if args['m2'] and not args['m1']:
		sys.exit("\nError: Must specify -1 and -2 if aligning paired end reads")
	# sanity check input values
	if args['mapid'] < 1 or args['mapid'] > 100:
		sys.exit("\nError: MAPQ must be between 1 and 100")
	if args['mapq'] < 0 or args['mapq'] > 100:
		sys.exit("\nError: MAPQ must be between 0 and 100")
	if args['baseq'] < 0 or args['baseq'] > 100:
		sys.exit("\nError: BASEQ must be between 0 and 100")

if __name__ == '__main__':

	program = get_program()
	args = get_arguments(program)
	check_arguments(program, args)
	create_directories(program, args)
	open_log(program, args)
	utility.print_copyright(args['log'])
	print_arguments(program, args)
	run_program(program, args)
