
def pileup_to_bases(pileup, ref_base):
	""" Takes pileup string from mpileup output. Returns the aligned bases """
	bases = ''
	index = 0
	while True:
		# parse pileup
		pileup_char = pileup[index]
		if pileup_char == '^': # skip '^' (start of a read segment) and ASCII character following '^'
			index += 2
		elif pileup_char == '$': # skip '$' (end of a read segment)
			index += 1
		elif pileup_char in ['+', '-']: # insertions and deletions (max 999 bp)
			for ndigits in [3,2,1]:
				indel_length = pileup[index+1:index+1+ndigits]
				if indel_length.isdigit(): break
			index += 1 + ndigits + int(indel_length)
		else: # keep positions for aligned bases
			if pileup_char in ['.',',']:
				bases += ref_base.upper()
			else:
				bases += pileup_char.upper()
			index += 1
		# stop when end of pileup reached
		if index == len(pileup):
			break
	return bases

def convert_from_ascii_quality(asciis):
	""" Convert quality scores to Sanger encoded (Phred+33) ascii values """
	ascii = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQR"""
	ascii_to_score = dict((y,x) for x,y in zip(range(0,50),list(ascii)))
	return [ascii_to_score[x] for x in asciis]

def parse_pileup(inpath):
	""" Yields formatted records from mpileup output file """
	for line in open(inpath):
		r = line.rstrip().split()
		ref_id = r[0]
		ref_pos = int(r[1])
		ref_base = r[2].upper()
		depth = int(r[3])
		if depth == 0:
			bases = ''
			phred_qualities = ''
			map_qualities = ''
		else:
			bases = pileup_to_bases(r[4], ref_base)
			phred_qualities = convert_from_ascii_quality(r[5])
			map_qualities = convert_from_ascii_quality(r[6])
		yield {'ref_id':ref_id, 'ref_pos':ref_pos,
		       'ref_base':ref_base, 'depth':depth,
			   'bases':bases, 'phred_qualities':phred_qualities,
			   'map_qualities':map_qualities}


inpath = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/example/ex1/pileup/58110.pileup'

for r in parse_pileup(inpath):
	print r











