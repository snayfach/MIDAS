#! /usr/bin/env python

# Parses pileup base string and returns the counts for all possible alleles

import csv, gzip, sys

class Pileup:
	def __init__(self, row):
		self.ref_id = row[0]
		self.ref_pos = int(row[1])
		self.ref_allele = row[2].upper()
		self.pileup = row[4]
		self.compute_counts()
		self.ref_count = self.counts[self.ref_allele]
		self.ref_freq = self.ref_count/float(self.depth) if self.depth > 0 else 0.0
		self.define_alt_allele()

	def compute_counts(self):
		""" get counts of 4 alleles """
		# For format specs, see: http://www.htslib.org/doc/samtools.html
		self.depth = 0 # count A,T,C,G,N
		self.index = 0 # pointer to current position in pileup string
		self.length = len(self.pileup)
		self.counts = {'A':0,'G':0,'C':0,'T':0,'N':0}
		while self.index < self.length:
			if self.pileup[self.index] in ['$', '*']:
				# $ denotes end of read segment, * denotes deletion
				self.index += 1
			elif self.pileup[self.index] == '^':
				# skip two characters when encountering '^' as it indicates
				#   a read start mark and the read mapping quality
				self.index += 2
			elif self.pileup[self.index] in ['+', '-']:
				# skip indels
				self.index += ( 2 + int(self.pileup[self.index+1]) )
			elif self.pileup[self.index] in ['.', ',']:
				# reference base
				self.counts[self.ref_allele] += 1
				self.depth += 1
				self.index += 1
			else:
				# non-reference base
				self.counts[self.pileup[self.index].upper()] += 1
				self.depth += 1
				self.index += 1

	def define_alt_allele(self):
		""" determine alternate allele """
		self.alt_allele = 'NA'
		self.alt_count = 0
		for allele in list('ATCG'):
			if allele == self.ref_allele:
				continue
			if self.counts[allele] > self.alt_count:
				self.alt_allele = allele
				self.alt_count = self.counts[allele]

	def allele_string(self):
		""" format allele counts """
		return ','.join([str(self.counts[_]) for _ in list('ATCG')])

def main(pileup_path):
	pileup_file = gzip.open(pileup_path)
	fnames =['ref_id', 'ref_pos', 'ref_allele', 'depth', 'pileup', 'qualities']
	csv.field_size_limit(sys.maxsize) # to deal with long fields
	pileup_reader = csv.reader(pileup_file, delimiter='\t')
	for row in pileup_reader:
		if row[2].upper() not in ['A','T','C','G']:
			continue
		else:
			yield Pileup(row)
	pileup_file.close()

if __name__ == '__main__':
    main(inpath)





