import os, sys, gzip, Bio.SeqIO, argparse, time, shutil
from multiprocessing import Process

def parse_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(usage='%s [options]' % os.path.basename(__file__))
	parser.add_argument('--indir', dest='indir', help='Input directory', required=True)
	parser.add_argument('--outdir', dest='outdir', help='Output directory', required=True)
	parser.add_argument('--genome_cluster', dest='genome_cluster', help='Genome-cluster ID', required=True)
	parser.add_argument('--sample_depth', dest='sample_depth', help='Minimum depth to include sample', type=float, default=10.0)
	parser.add_argument('--snp_depth', dest='snp_depth', help='Minimum relative-depth to include SNP', type=float, default=0.30)
	parser.add_argument('--threads', dest='threads', help='Number of threads to use', type=int, default=1)
	return vars(parser.parse_args())

def get_marker_cov(inpath):
	""" Get marker gene coverage from pangenome coverage file """
	for line in gzip.open(inpath):
		cov, copy = [float(x) for x in line.rstrip().split()[1:3]]
		if cov > 0 and copy > 0:
			return cov/copy
	return 0

def parse_snps(inpath):
	""" Yield formatted record from snps file """
	f = [('ref_id',str), ('ref_pos',str),
		 ('ref_allele',str), ('alt_allele',str),
		 ('cons_allele',str), ('count_alleles',int),
		 ('count_ref',int), ('count_alt',int),
		 ('depth',int), ('ref_freq',float)]
	infile = gzip.open(inpath)
	next(infile)
	for line in infile:
		r = line.rstrip().split()
		if int(r[-2]) == 0: continue
		yield dict([(f[i][0], f[i][1](j)) for i,j in enumerate(r)])

def select_sample(sample_id, args):
	""" Determine whether sample has sufficient depth for target genome-cluster """
	snps_path = '/'.join([args['indir'], sample_id, 'snps', '%s.snps.gz' % args['genome_cluster']])
	cov_path = '/'.join([args['indir'], sample_id, 'coverage', '%s.cov.gz' % args['genome_cluster']])
	if not all([os.path.isfile(x) for x in [snps_path, cov_path]]): # input file(s) doesn't exist
		return False
	elif get_marker_cov(cov_path) < args['sample_depth']: # insufficient depth
		return False
	else:
		return True

def identify_samples(args):
	""" Identify samples with sufficient depth for target genome-cluster """
	samples = []
	for sample_id in os.listdir(args['indir']):
		if select_sample(sample_id, args):
			samples.append(sample_id)
	return samples

def write_snp_list(args, sample_id):
	""" Extract list of SNPs with sufficient depth for target sample_id, genome_cluster """
	# output paths
	outdir = '/'.join([args['outdir'], 'snp_lists'])
	out_file = gzip.open('/'.join([outdir, '%s.snp_list.gz' % sample_id]), 'w')
	# input paths
	snps_path = '/'.join([args['indir'], sample_id, 'snps', '%s.snps.gz' % args['genome_cluster']])
	cov_path = '/'.join([args['indir'], sample_id, 'coverage', '%s.cov.gz' % args['genome_cluster']])
	# id and write snps to out_file
	marker_cov = get_marker_cov(cov_path)
	infile = gzip.open(snps_path)
	next(infile)
	for line in infile:
		r = line.split()
		depth = float(r[-2])
		if depth/marker_cov >= args['snp_depth']:
			out_rec = [str(x) for x in [r[0], r[1]]]
			out_file.write('\t'.join(out_rec)+'\n')

def parallel_process(function, args_list, threads):
	""" Run function using multiple threads """
	processes = []
	for pargs in args_list: # run function for each set of args in args_list
		p = Process(target=function, kwargs=pargs)
		processes.append(p)
		p.start()
		# control number of active processes
		while len(processes) >= threads:
			time.sleep(1)
			# remove processes that are no longer alive
			indexes = []
			for index, process in enumerate(processes):
				if process.is_alive(): indexes.append(index)
			processes = [processes[i] for i in indexes]
	# wait until there are no active processes
	while len(processes) > 0:
		time.sleep(1)
		indexes = []
		# remove processes that are no longer alive
		for index, process in enumerate(processes):
			if process.is_alive(): indexes.append(index)
		processes = [processes[i] for i in indexes]

def merge_snps(args, samples, min_prev):
	""" Select SNPs that meet/exceed minimum prevalence """
	# count the number of times each SNP observed across samples
	snp_counts = {}
	indir = '/'.join([args['outdir'], 'snp_lists'])
	for sample_id in samples:
		infile = gzip.open('/'.join([indir, '%s.snp_list.gz' % sample_id]))
		for line in infile:
			chr, pos = line.rstrip().split()
			snp = (chr, pos)
			if snp in snp_counts:
				snp_counts[snp] += 1
			else:
				snp_counts[snp] = 1
	# select SNPs that meet/exceed minimum prevalence
	merged_snps = set([])
	count_samples = len(samples)
	for snp, count in snp_counts.items():
		if float(count)/count_samples >= min_prev:
			merged_snps.add(snp)
	return merged_snps

def write_results(args, samples, snps):
	""" Write consensus sequences, SNP depth, and SNP frequencies to disk for specified samples and SNPs """
	# Open output files
	cons_file = gzip.open('/'.join([args['outdir'], '%s.consensus.gz' % args['genome_cluster']]), 'w')
	cov_file = gzip.open('/'.join([args['outdir'], '%s.coverage.gz' % args['genome_cluster']]), 'w')
	freq_file = gzip.open('/'.join([args['outdir'], '%s.ref_freq.gz' % args['genome_cluster']]), 'w')
	# Write data
	for sample_id in samples:
		cons_file.write('>'+sample_id+'\n')
		cov_file.write(sample_id)
		freq_file.write(sample_id)
		inpath = '/'.join([args['indir'], sample_id, 'snps', '%s.snps.gz' % args['genome_cluster']])
		for r in parse_snps(inpath):
			snp = (r['ref_id'], r['ref_pos'])
			if snp in snps:
				cons_file.write(r['cons_allele'])
				cov_file.write('\t%s' % r['depth'])
				freq_file.write('\t%s' % r['ref_freq'])
		cons_file.write('\n')
		cov_file.write('\n')
		freq_file.write('\n')

if __name__ == '__main__':

	args = parse_arguments()
	
	# create temp dirs
	tmpdir = '/'.join([args['outdir'], 'snp_lists'])
	if not os.path.isdir(args['outdir']): os.mkdir(args['outdir'])
	if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
	
	# id samples with sufficient depth
	print("Identifying samples")
	samples = identify_samples(args) # ['SRR1761671', 'SRR1747043', 'SRR1761673', 'SRR1747031'] # 

	# id snps with sufficient depth
	print("Identifying SNPs")
	parallel_process(write_snp_list, [{'args':args,'sample_id':x} for x in samples], args['threads'])

	# merge snps
	print("Intersecting SNPs")
	snps = merge_snps(args, samples, min_prev=1.0)

	# write results
	print("Writing results for intersected SNPs")
	write_results(args, samples, snps)

	# cleanup
	shutil.rmtree(tmpdir)

