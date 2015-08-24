import pysam, sys

def compute_perc_id(aln):
	""" Compute percent identity for read """
	length = aln.query_length
	edit = dict(aln.tags)['NM']
	return 100 * (length - edit)/float(length)

def filter_bam(inpath, outpath, pid):
	""" Filter records from bamfile and write to temporary output file """
	infile = pysam.AlignmentFile(inpath, 'rb')
	outfile = pysam.AlignmentFile(outpath, 'wb', template=infile)
	for aln in infile:
		if compute_perc_id(aln) >= pid:
			outfile.write(aln)

if __name__ == '__main__':
	filter_bam(sys.argv[1], sys.argv[2], float(sys.argv[3]))