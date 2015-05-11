""" Map pangene ids to phyeco ids in order to normalize pangene coverage estimates """

import os, gzip

project_dir = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation'

gene_to_phyeco = {}
inpath = os.path.join(project_dir, 'genome_clustering2/db/marker_to_gene.txt')
for line in open(inpath):
	phyeco_id, gene_id = line.rstrip().split()
	gene_to_phyeco[gene_id] = phyeco_id

clusters_dir = os.path.join(project_dir, 'cnv_detection2/microbe_cnv/genome_clusters')

for cluster in os.listdir(clusters_dir):
	outfile = gzip.open('/'.join([clusters_dir, cluster, 'pangene_to_phyeco.gz' ]), 'w')
	outfile.write('\t'.join(['pangene_id', 'type', 'phyeco_id'])+'\n')
	bedfile = open('/'.join([clusters_dir, cluster, '%s.bed' % cluster]))
	pangene_to_phyeco = set([])
	for line in bedfile:
		sid, start, end, gene_id, pangene_id = line.rstrip().split()
		if gene_id in gene_to_phyeco:
			pangene_to_phyeco.add((pangene_id, gene_to_phyeco[gene_id]))
	for record in pangene_to_phyeco:
		outfile.write('\t'.join(record)+'\n')


