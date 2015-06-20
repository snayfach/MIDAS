import os, sys, gzip, Bio.SeqIO

indir = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/microbe_cnv_out'
outdir = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/snp_profiles'
my_cluster_id = sys.argv[1] # = ["62236", "57453", "61481"]

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

# Aggregate run_ids and file paths across clusters
print("Aggregating runs by genome-cluster")
cluster_to_runs = {}
for run in os.listdir(indir):
	if not os.path.isdir('/'.join([indir, run, 'snps'])):
		print ("\twarning, no snps dir for %s" % run)
		continue
	for file in os.listdir('/'.join([indir, run, 'snps'])):
		cluster_id = file.split('.')[0]
		if cluster_id != my_cluster_id: continue
		marker_cov = get_marker_cov('/'.join([indir, run, 'coverage', '%s.cov.gz' % cluster_id]))
		if marker_cov < 3.0:
			continue
		elif cluster_id not in cluster_to_runs:
			cluster_to_runs[cluster_id] = [[run, marker_cov]]
		else:
			cluster_to_runs[cluster_id].append([run, marker_cov])

# Build snp matrices
for cluster_id, runs in cluster_to_runs.items():

	# skip off-target clusters
	if cluster_id != my_cluster_id: continue

	# store reference positions
	ref_alleles = {}

	# store data across samples
	print("Storing snps for %s runs from genome-cluster %s" % (len(runs), cluster_id))
	run_to_snps = {}
	for run, marker_cov in runs:
		print("\t%s" % run)
		run_to_snps[run] = {} # init snps for sample
		snps_path = '/'.join([indir, run, 'snps', '%s.snps.gz' % cluster_id])
		for r in parse_snps(snps_path):
			pos = (r['ref_id'], int(r['ref_pos']))
			if r['depth']/marker_cov >= 1/3:
				run_to_snps[run][pos] = {'base':r['cons_allele'], 'cov':r['depth'], 'freq':r['ref_freq']}
				ref_alleles[pos] = r['ref_allele']

	# merge snps
	print("Intersecting snps")
	positions = set(run_to_snps.values()[0])
	for snps in run_to_snps.values()[1:]:
		positions = positions.intersection(set(snps))
	positions = sorted(positions)
	print("\t%s positions")
	
	# open output files
	cons_file = gzip.open(os.path.join(outdir, '%s.consensus.gz' % cluster_id), 'w')
	cov_file = gzip.open(os.path.join(outdir, '%s.coverage.gz' % cluster_id), 'w')
	freq_file = gzip.open(os.path.join(outdir, '%s.ref_freq.gz' % cluster_id), 'w')
	
	# write headers
	cov_file.write('run'); freq_file.write('run')
	for pos in positions:
		cov_file.write('\t'); cov_file.write('|'.join([str(x) for x in pos]))
		freq_file.write('\t'); freq_file.write('|'.join([str(x) for x in pos]))
	cov_file.write('\n'); freq_file.write('\n')

	# write reference sequence
	cons_file.write('>'+'reference'+'\n')
	for pos in positions:
		cons_file.write(ref_alleles[pos])
	cons_file.write('\n')

	# write results
	print("Writing intersected snps")
	for run in run_to_snps:
		cons_file.write('>'+run+'\n'); cov_file.write(run); freq_file.write(run)
		for i, pos in enumerate(positions):
			snp = run_to_snps[run][pos]
			cons_file.write(snp['base'])
			cov_file.write('\t'); cov_file.write(str(snp['cov']))
			freq_file.write('\t'); freq_file.write(str(snp['freq']))
		cons_file.write('\n'); cov_file.write('\n'); freq_file.write('\n')

