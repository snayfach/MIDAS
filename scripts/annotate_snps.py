#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

## TO DO
#  -Handle rna genes

__version__ = '0.0.2'

import argparse, sys, os, gzip, Bio.SeqIO

def parse_arguments():
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
	parser.add_argument('-i', '--snp_list', dest='in', type=str, help='SNP list', required=True)
	parser.add_argument('-o', '--annotations', dest='out', type=str, help='SNP annotations', required=True)
	parser.add_argument('-m', '--max_snps', dest='max_snps', type=int, help='Maximum # of SNPs to annotate (all)')
	parser.add_argument('-v', '--verbose', action='store_true', default=False)
	return vars(parser.parse_args())


def parse_snps(inpath):
	infile = open(inpath)
	fields = next(infile).rstrip().split()
	for line in infile:
		values = line.rstrip().split()
		record = dict([(f,v) for f,v in zip(fields, values)])
		record['ref_pos'] = int(record['ref_pos'])
		yield record

def read_genome():
	cluster_id = os.path.basename(args['in']).split('.')[0]
	pkg_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	inpath = '%s/ref_db/genome_clusters/%s/representative.fna.gz' % (pkg_dir, cluster_id)
	infile = gzip.open(inpath)
	genome = {}
	for r in Bio.SeqIO.parse(infile, 'fasta'):
		genome[r.id] = r.seq
	return genome

def read_genes(skip_rna=True):
	""" Read in gene coordinates from features file """
	cluster_id = os.path.basename(args['in']).split('.')[0]
	pkg_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	genes_path = '%s/ref_db/genome_clusters/%s/representative.features.gz' % (pkg_dir, cluster_id)
	genes = []
	infile = gzip.open(genes_path)
	fields = next(infile).rstrip().split()
	for line in infile:
		values = line.rstrip().split()
		gene = dict([(f,v) for f,v in zip(fields, values)])
		if skip_rna and gene['gene_type'] == 'rna':
			continue
		else:
			gene['accession'] = 'accn|%s' % gene['accession']
			gene['start'] = int(gene['start'])
			gene['end'] = int(gene['end'])
			genes.append(gene)
	return genes

def check_snps():
	last_snp = {'ref_id':'accn|NZ_DS499510'}
	for snp in parse_snps(args['in']):
		if snp['ref_id'] < last_snp['ref_id']:
			sys.exit("Accessions not sorted")

def check_genes(genes, genomes):
	last_gene = {'accession':'accn|NZ_DS499510'}
	for gene in genes:
		if not gene['accession'] in genome:
			sys.exit("Incorrect scaffold id")
		if (gene['end'] - gene['start'] + 1) % 3 != 0:
			sys.exit("Incorrect gene length")
		if gene['accession'] < last_gene['accession']:
			sys.exit("Accessions not sorted")
		last_gene = gene

def rev_comp(seq):
	""" Reverse complement sequence """
	d = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
	return(''.join([d[i] for i in list(seq[::-1])]))

def get_gene_seq(gene, genome):
	""" Fetch nucleotide sequence of gene from genome """
	seq = genome[gene['accession']][gene['start']-1:gene['end']].upper()
	if gene['strand'] == '-':
		return(rev_comp(seq))
	else:
		return(seq)

def translate(codon):
	""" Translate individual codon """
	codontable = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
	'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
	}
	return codontable[str(codon)]

def index_replace(x, y, i):
	""" Replace character at index i  in string x with y"""
	z = list(x)
	z[i] = y
	return(''.join(z))

def classify_site(ref_codon, codon_pos):
	""" Classify coding site as ND, 2D, 3D, or 4D """
	count = 0
	degeneracy = {0:'ND',1:'2D',2:'3D',3:'4D'}
	ref_aa = translate(ref_codon)
	ref_allele = list(ref_codon)[codon_pos]
	for allele in ['A','T','C','G']:
		if allele == ref_allele:
			continue
		elif translate(index_replace(ref_codon, allele, codon_pos)) == ref_aa:
			count += 1
	return degeneracy[count]

def classify_snp(ref_codon, alt_allele, codon_pos):
	""" Classify SNP an SYN or NS """
	alt_codon = index_replace(ref_codon, alt_allele, codon_pos)
	alt_aa = translate(alt_codon)
	if translate(ref_codon) == alt_aa:
		return 'SYN'
	else:
		return 'NS'

def annotate_site_and_snp(snp, genome):
	""" Annotate variant and reference site """
	global genes
	site_type = None # NC, NS, 2D, 3D, 4D
	snp_type = None # SYN, NS, NA
	gene = None
	while True:
		gene = genes[0]
		# snp downstream of last gene
		if len(genes) == 0:
			site_type = 'NC'
			snp_type = 'NA'
			gene = 'NA'
			break
		# snp upstream of next gene
		elif (snp['ref_id'] < gene['accession'] or
			 (snp['ref_id'] == gene['accession'] and
			  snp['ref_pos'] < gene['start'])):
			site_type = 'NC'
			snp_type = 'NA'
			gene = 'NA'
			break
		# snp downstream previous gene
		elif (snp['ref_id'] > gene['accession'] or
			  (snp['ref_id'] == gene['accession'] and
			  snp['ref_pos'] > gene['end'])):
			genes = genes[1:]
		# snp in gene
		else:
			gene_pos = snp['ref_pos']-gene['start'] if gene['strand']=='+' else gene['end']-snp['ref_pos'] # position of snp in gene
			codon_pos=gene_pos%3 # position of snp in codon
			seq = get_gene_seq(gene, genome) # gene sequence (oriented start to stop)
			ref_codon = [seq[i:i+3] for i in range(0, len(seq), 3)][gene_pos/3]
			site_type = classify_site(ref_codon, codon_pos)
			if snp['alt_allele'] == 'NA': # fixed reference allele
				snp_type = 'NA'
			else:
				alt_allele = snp['alt_allele'] if gene['strand'] == '+' else rev_comp(snp['alt_allele'])
				snp_type = classify_snp(ref_codon, alt_allele, codon_pos)
			break
	# record annotations
	snp['site_type'] = site_type
	snp['snp_type'] = snp_type
	snp['gene_id'] = gene['gene_id'].split('|')[-1] if gene != 'NA' else 'NA'


def write_annotation(snp, outfile, fields):
	""" Write record to output file """
	values = [str(x) for x in [snp[field] for field in fields]]
	outfile.write('\t'.join(values)+'\n')

if __name__ == '__main__':
	args = parse_arguments()
	genome = read_genome()
	genes = read_genes()
	outfile = open(args['out'], 'w')
	fields = ['ref_id', 'ref_pos', 'count_alt', 'alt_allele', 'site_type', 'snp_type', 'gene_id']
	outfile.write('\t'.join(fields)+'\n')
	for i, snp in enumerate(parse_snps(args['in'])):
		annotate_site_and_snp(snp, genome)
		write_annotation(snp, outfile, fields)
		if args['max_snps'] and i >= args['max_snps']: break

	




