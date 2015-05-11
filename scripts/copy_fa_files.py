# This script will copy fasta files for representative reference genomes for each genome-cluster
# Also, genomes will be indexed using samtools faidx

# libs
import os, gzip, shutil, subprocess

# dirs
work_dir = '/mnt/data/work/pollardlab/snayfach'
genomes_dir = '/'.join([work_dir, 'genomes/PATRIC3/ftp.patricbrc.org/patric2/patric3/genomes'])
gc_dir = '/'.join([work_dir, 'projects/strain_variation/microbe_cnv/genome_clusters'])

samtools = '/'.join([work_dir, 'projects/strain_variation/microbe_cnv/lib/samtools-1.1/samtools'])

# read in representative genomes
gc_to_rep = {}
infile = open('/'.join([work_dir, 'projects/strain_variation/genome_clustering2/clusters_combined/top_30/0.035/cluster_to_centroid.txt']))
for line in infile:
	cluster_id, genome_id = line.rstrip().split()
	gc_to_rep[cluster_id] = genome_id

# loop over genome-cluster directories
for index, gc in enumerate(os.listdir(gc_dir)):
	
	print index
	centroid_id = gc_to_rep[gc]
	
	# copy representative genome
	src = '/'.join([genomes_dir, centroid_id, '%s.fna' % centroid_id])
	dest = '/'.join([gc_dir, gc, 'representative.fna'])
	shutil.copyfile(src, dest)

	# index fasta file
	process = subprocess.Popen("%s faidx %s" % (samtools, dest), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()



