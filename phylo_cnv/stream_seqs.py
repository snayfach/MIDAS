#!/usr/bin/env python

import os, sys, gzip, argparse

def iopen(inpath):
	""" Open input file for reading regardless of compression [gzip, bzip] or python version """
	ext = inpath.split('.')[-1]
	# Python2
	if sys.version_info[0] == 2:
		if ext == 'gz': return gzip.open(inpath)
		elif ext == 'bz2': return bz2.BZ2File(inpath)
		else: return open(inpath)
	# Python3
	elif sys.version_info[0] == 3:
		if ext == 'gz': return io.TextIOWrapper(gzip.open(inpath))
		elif ext == 'bz2': return bz2.BZ2File(inpath)
		else: return open(inpath)

def readfq(fp):
	""" https://github.com/lh3/readfq/blob/master/readfq.py 
		A generator function for parsing fasta/fastq records """
	last = None # this is a buffer keeping the last unprocessed line
	while True: # mimic closure; is it a bad idea?
		if not last: # the first record or a record following a fastq
			for l in fp: # search for the start of the next record
				if l[0] in '>@': # fasta/q header line
					last = l[:-1] # save this line
					break
		if not last: break
		name, seqs, last = last[1:].partition(" ")[0], [], None
		for l in fp: # read the sequence
			if l[0] in '@+>':
				last = l[:-1]
				break
			seqs.append(l[:-1])
		if not last or last[0] != '+': # this is a fasta record
			yield name, ''.join(seqs), None # yield a fasta record
			if not last: break
		else: # this is a fastq record
			seq, leng, seqs = ''.join(seqs), 0, []
			for l in fp: # read the quality
				seqs.append(l[:-1])
				leng += len(l) - 1
				if leng >= len(seq): # have read enough quality
					last = None
					yield name, seq, ''.join(seqs); # yield a fastq record
					break
			if last: # reach EOF before reading enough quality
				yield name, seq, None # yield a fasta record instead
				break

def main():
	""" Run main pipeline """
	args = parse_args()
	reads = 0
	bp = 0
	for inpath in args['input']:
		infile = iopen(inpath)
		for name, seq, qual in readfq(infile):
			seq_len = len(seq)
			if args['read_length']: # trim/filter reads
				if seq_len < args['read_length']:
					continue
				else:
					seq = seq[0:args['read_length']]
					seq_len = len(seq)
			sys.stdout.write('>%s_%s\n%s\n' % (name, seq_len, seq))
			reads += 1
			bp += seq_len
			if reads == args['max_reads']:
				sys.stderr.write('%s\t%s' % (reads, bp)) # write number of reads, bp to stderr
				return
	sys.stderr.write('%s\t%s' % (reads, bp)) # write number of reads, bp to stderr

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-1', type=str, dest='m1')
	parser.add_argument('-2', type=str, dest='m2')
	parser.add_argument('-l', type=int, dest='read_length')
	parser.add_argument('-n', type=int, dest='max_reads', default=float('Inf'))
	args = vars(parser.parse_args())
	args['input'] = [args['m1']]
	if args['m2']: args['input'].append(args['m2'])
	return args

if __name__ == "__main__":
	main()



