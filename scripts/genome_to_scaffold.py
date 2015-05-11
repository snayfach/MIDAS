# libs
import os, gzip

# dirs
work_dir = '/mnt/data/work/pollardlab/snayfach'
genomes_dir = '/'.join([work_dir, 'genomes/PATRIC3/ftp.patricbrc.org/patric2/patric3/genomes'])
gc_dir = '/'.join([work_dir, 'projects/strain_variation/cnv_detection2/microbe_cnv/genome_clusters'])

# read in map files
gc_to_genomes = {}
infile = open('/'.join([work_dir, 'projects/strain_variation/genome_clustering2/clusters_combined/top_30/0.035/genome_clusters.txt']))
for line in infile:
	cluster_id, genomes = line.rstrip().split()
	gc_to_genomes[cluster_id] = genomes.split(',')

# map genome id to list of scaffold ids
for gc in os.listdir(gc_dir):
	outfile = gzip.open('/'.join([gc_dir, gc, 'genome_to_scaffold.gz']), 'w')
	for genome_id in gc_to_genomes[gc]:
		inpath = '/'.join([genomes_dir, genome_id, '%s.fna' % genome_id])
		for line in open(inpath):
			if line[0] != '>': continue
			sid = line.rstrip().lstrip('>').split()[0].split('|')[1]
			outfile.write('\t'.join([genome_id, sid])+'\n')