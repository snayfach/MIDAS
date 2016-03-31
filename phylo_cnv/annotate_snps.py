#!/usr/bin/env python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

## TO DO
#  -Handle rna genes

import argparse, sys, os, gzip, Bio.SeqIO
from phylo_cnv import utility

def read_genome(db, species_id):
	""" Read in representative genome from reference database """
	inpath = '%s/genome_clusters/%s/representative.fna.gz' % (db, species_id)
	infile = gzip.open(inpath)
	genome = {}
	for r in Bio.SeqIO.parse(infile, 'fasta'):
		genome[r.id] = r.seq
	infile.close()
	return genome

def read_genes(db, species_id, skip_rna=True):
	""" Read in gene coordinates from features file """
	genes_path = '%s/genome_clusters/%s/representative.features.gz' % (db, species_id)
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
	infile.close()
	return genes

def check_snps():
	""" Check that accessions are sorted """
	last_snp = {'ref_id':'accn|NZ_DS499510'}
	for snp in utility.parse_file(args['in']):
		if snp['ref_id'] < last_snp['ref_id']:
			sys.exit("Accessions not sorted")

def check_genes(genes, genomes):
	""" Check that accessions are sorted and gene lengths are correct """
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
	""" Replace character at index i in string x with y"""
	z = list(x)
	z[i] = y
	return(''.join(z))

def classify_site(ref_codon, codon_pos):
	""" Classify coding site as ND, 2D, 3D, or 4D """
	count = 0
	degeneracy = {0:'1D',1:'2D',2:'3D',3:'4D'}
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

def annotate_site(site, genes, genome):
	""" Annotate variant and reference site """
	site['ref_pos'] = int(site['ref_pos'])
	site['gene_id'] = 'NA'
	site['snp_type'] = {'A':'NA','T':'NA','C':'NA','G':'NA'}
	while True:
		# snp downstream of last gene
		if len(genes) == 0:
			site['site_type'] = 'NC'
			gene = 'NA'
			break
		else:
			gene = genes[0]
		# snp upstream of next gene
		if (site['ref_id'] < gene['accession'] or
			 (site['ref_id'] == gene['accession'] and
			  site['ref_pos'] < gene['start'])):
			site['site_type'] = 'NC'
			gene = 'NA'
			break
		# snp downstream previous gene
		elif (site['ref_id'] > gene['accession'] or
			  (site['ref_id'] == gene['accession'] and
			  site['ref_pos'] > gene['end'])):
			genes = genes[1:]
		# snp in gene
		else:
			ref_codon, codon_pos = fetch_ref_codon(site, gene, genome)
			site['gene_id'] = gene['gene_id'].split('|')[-1]
			site['site_type'] = classify_site(ref_codon, codon_pos)
			for alt_allele in ['A','T','C','G']:
				site['snp_type'][alt_allele] = classify_snp(ref_codon, alt_allele, codon_pos)
			break

def fetch_ref_codon(site, gene, genome):
	""" Fetch codon within gene for given site """
	gene_pos = site['ref_pos']-gene['start'] if gene['strand']=='+' else gene['end']-site['ref_pos'] # position of snp in gene
	codon_pos=gene_pos%3 # position of snp in codon
	seq = get_gene_seq(gene, genome) # gene sequence (oriented start to stop)
	ref_codon = [seq[i:i+3] for i in range(0, len(seq), 3)][gene_pos/3]
	return ref_codon, codon_pos

def write_record(site, outfile):
	""" Write site info to output file """
	values = []
	for field in ['ref_id', 'ref_pos', 'ref_allele', 'gene_id', 'site_type']:
		values.append(str(site[field]))
	for snp in ['A','T','C','G']:
		values.append(site['snp_type'][snp])
	outfile.write('\t'.join(values)+'\n')

def open_outfile(args):
	""" Open snp_info file for writing """
	outfile = open('%s/%s.snps.info' % (args['outdir'], args['species_id']), 'w')
	fields = ['ref_id', 'ref_pos', 'ref_allele', 'gene_id', 'site_type', 'snp_A', 'snp_T', 'snp_C', 'snp_G']
	outfile.write('\t'.join(fields)+'\n')
	return outfile

def main(args):
	genome = read_genome(args['db'], args['species_id'])
	genes = read_genes(args['db'], args['species_id'])
	snpinfo = open_outfile(args)
	snplist = '%s/%s.snps.list' % (args['outdir'], args['species_id'])
	for i, site in enumerate(utility.parse_file(snplist)):
		annotate_site(site, genes, genome)
		write_record(site, snpinfo)
		if args['max_sites'] and i >= args['max_sites']: break
	snpinfo.close()

	




