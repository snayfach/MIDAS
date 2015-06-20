import os, sys, gzip

indir = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/microbe_cnv_out'
outdir = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation/microbe_cnv/application/pangenome_profiles'

def get_marker_cov(inpath):
	""" Get marker gene coverage from pangenome coverage file """
	for line in gzip.open(inpath):
		cov, copy = [float(x) for x in line.rstrip().split()[1:3]]
		if cov > 0 and copy > 0:
			return cov/copy
	return 0

# Aggregate run_ids, marker coverages, and file paths across clusters
cluster_to_paths = {}
for run in os.listdir(indir):
	cov_dir = '/'.join([indir, run, 'coverage'])
	if not os.path.isdir(cov_dir):
		print ("Warning no coverage dir for %s" % run)
		continue
	for file in os.listdir(cov_dir):
		cluster_id = file.split('.')[0]
		covpath = '/'.join([cov_dir, file])
		if cluster_id not in cluster_to_paths:
			cluster_to_paths[cluster_id] = [(run, covpath)]
		else:
			cluster_to_paths[cluster_id].append((run, covpath))

# Build pangenome matrices
for cluster_id, paths in cluster_to_paths.items():

	# open out for cluster_id
	outpath = os.path.join(outdir, '%s.pg_profile.gz' % cluster_id)
	outfile = gzip.open(outpath, 'w')
	print ("Now writing %s" % cluster_id)
	
	# write header
	pangene_ids = [_.rstrip().split()[0] for _ in gzip.open(paths[0][1])]
	field_names = ['run', 'marker_coverage'] + pangene_ids
	outfile.write('\t'.join(field_names)+'\n')
	
	# write rest of lines
	for run, path in paths:
		marker_cov = get_marker_cov(path)
		outfile.write(run+'\t'+str(marker_cov))
		for line in gzip.open(path):
			copy_number = line.rstrip().split()[-1]
			outfile.write('\t')
			outfile.write(copy_number)
		outfile.write('\n')
	outfile.close()


