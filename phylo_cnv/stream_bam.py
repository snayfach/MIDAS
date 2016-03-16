import pysam, sys, numpy as np

def compute_perc_id(aln):
	""" Compute percent identity of aligned region on read """
	length = len(aln.query_alignment_sequence)
	edit = dict(aln.tags)['NM']
	return 100 * (length - edit)/float(length)

def compute_aln_cov(aln):
	""" Compute percent identity for paired-end read """
	aln_cov = len(aln.query_alignment_sequence)/float(aln.query_length)
	return aln_cov

def compute_read_qual(aln):
	""" average read qualiy """
	np.mean(aln.query_qualities)
	return aln_cov

def find_indexes(string):
	""" Find indexes of all N's in string """
	return [i for i, s in enumerate(string) if s == 'N']

def fix_qualities(aln, qual=2):
	""" Fix base quality of N-calls """
	indexes = find_indexes(aln.query_sequence)
	if len(indexes) == 0:
		return
	else:
		qualities = aln.query_qualities
		for i in indexes:
			qualities[i] = qual
		aln.query_qualities = qualities

def filter_bam(inpath, outpath, pid, min_qual):
	""" Filter records from bamfile and write to temporary output file """
	infile = pysam.AlignmentFile(inpath, 'rb')
	outfile = pysam.AlignmentFile(outpath, 'wb', template=infile)
	for aln in infile:
		fix_qualities(aln)
		if compute_perc_id(aln) < pid:
			continue
		elif np.mean(aln.query_qualities) < min_qual:
			continue
		else:
			outfile.write(aln)

if __name__ == '__main__':

	filter_bam(inpath=sys.argv[1], outpath=sys.argv[2], pid=float(sys.argv[3]), min_qual=float(sys.argv[4]))





