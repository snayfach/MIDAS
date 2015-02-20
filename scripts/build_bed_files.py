import os, gzip

# in/out directories
bedfiles_dir = '/mnt/data/work/pollardlab/snayfach/genomes/PATRIC/bed_files'
pangenomes_dir = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation/pan_genomics/pan_genomes/0.90'
out_dir = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation/cnv_detection/develop_software/genome_clusters'

i = 0

for cluster_id in os.listdir(pangenomes_dir):

	# open outfile
	outfile = gzip.open('/'.join([out_dir, cluster_id, cluster_id+'.bed.gz']), 'w')
	header = ['genome_id', 'pangene_id', 'type', 'gene_id', 'scaffold_id', 'start', 'end']
	outfile.write('\t'.join(header)+'\n')

	# gene to pangene
	gene_to_pgene = {}
	infile = gzip.open('/'.join([pangenomes_dir, cluster_id, 'gene_to_cluster.txt.gz']))
	next(infile)
	for line in infile:
		pgene_id, type, gene_id, genome_id = line.rstrip().split()
		gene_to_pgene[gene_id] = (pgene_id, type)
	
	# gene to genomic location(s)
	infile = gzip.open('/'.join([pangenomes_dir, cluster_id, 'genomes.txt.gz']))
	for line in infile:
		i += 1
		if not i % 100 : print i
		genome_id = line.rstrip()
		bedfile = gzip.open(os.path.join(bedfiles_dir, genome_id+'.bed.gz'))
		for line in bedfile:
			try:
				target_id, start, end, gene_id = line.rstrip().split()
				pgene_id, type = gene_to_pgene[gene_id]
				record = [genome_id, pgene_id, type, gene_id, target_id, start, end]
				outfile.write('\t'.join(record)+'\n')
			except Exception:
				print '\t', cluster_id, genome_id, gene_id
				pass