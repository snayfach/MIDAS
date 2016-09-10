#! /usr/bin/env python

# Parses pileup base string and returns the counts for all possible alleles

def parse_pileup(ref_allele, pileup):
	counts = {'A':0,'G':0,'C':0,'T':0,'-':[],'+':[]}
	pileup = pileup.upper()
	pileup = pileup.replace('$', '')
	pileup = pileup.replace('*', '')
	while pileup != '':
		if pileup[0] == '^':
			# skip two characters when encountering '^' as it indicates
			# a read start mark and the read mapping quality
			pileup = pileup[2:]
		elif pileup[0] in ['+', '-']:
			# skip indels
			pileup = pileup[2+int(pileup[1]):]
		elif pileup[0] in ['.', ',']:
			# reference base
			counts[ref_allele] += 1
			pileup = pileup[1:]
		else:
			# non-reference base
			counts[pileup[0]] += 1
			pileup = pileup[1:]
	return counts

def define_alt(ref_allele, counts):
	alt_allele = 'NA'
	alt_count = 0
	for allele in list('ATCG'):
		if allele == ref_allele:
			continue
		if counts[allele] > alt_count:
			alt_allele = allele
			alt_count = counts[allele]
	return alt_allele

def count_alleles(counts):
	return ','.join([str(counts[_]) for _ in list('ATCG')])

def ref_freq(ref_allele, counts, depth):
	return counts[ref_allele]/float(depth) if depth > 0 else 0.0

def main(infile):
	for line in infile:
		v = line.strip('\n').split('\t')
		r = {}
		r['ref_id'] = v[0]
		r['ref_pos'] = v[1]
		r['ref_allele'] = v[2].upper()
		r['depth'] = int(v[3])
		r['pileup'] = v[4]
		r['counts'] = parse_pileup(r['ref_allele'], r['pileup'])
		r['count_atcg'] = count_alleles(r['counts'])
		r['ref_freq'] = ref_freq(r['ref_allele'], r['counts'], r['depth'])
		r['alt_allele'] = define_alt(r['ref_allele'], r['counts'])
		yield r
	infile.close()

if __name__ == '__main__':
    main(inpath)





