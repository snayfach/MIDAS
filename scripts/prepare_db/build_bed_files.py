# import libs
import os, gzip

# in/out directories
project_dir = '/mnt/data/work/pollardlab/snayfach/projects/strain_variation'
genomes_dir = '/mnt/data/work/pollardlab/snayfach/genomes/PATRIC3/ftp.patricbrc.org/patric2/patric3/genomes'
out_dir = '/'.join([project_dir, 'cnv_detection2', 'microbe_cnv', 'genome_clusters'])

# loop over genome-clusters
cluster_ids = os.listdir('/'.join([project_dir, 'pan_genomics2', 'pan_genomes']))
for index, cluster_id in enumerate(cluster_ids):
	print index
	
	# open output file
	outfile = open('/'.join([out_dir, cluster_id, '%s.bed' % cluster_id]), 'w')
	
	# read in genomes
	inpath = '/'.join([project_dir, 'pan_genomics2', 'nr_genomes', cluster_id, 'representatives.txt'])
	if not os.path.isfile(inpath):
		print("%s not found" % os.path.basename(inpath) )
		continue
	genomes = [line.rstrip().split()[1] for line in open(inpath)]
	
	# map gene to pangene id
	gene_to_pangene = {}
	inpath = '/'.join([project_dir, 'pan_genomics2', 'pan_genomes', cluster_id, 'gene_to_cluster.txt'])
	if not os.path.isfile(inpath):
		print("%s not found" % os.path.basename(inpath) )
		continue
	infile = open(inpath)
	next(infile)
	for line in infile:
		clust_id, clust_type, gene_id, genome_id = line.rstrip().split()
		pangene_id = '_'.join([clust_id, clust_type])
		gene_to_pangene[gene_id] = pangene_id

	# loop over genomes
	for genome_id in genomes:
	
		for type in ['cds', 'rna']:
		
			# check feature file exists
			inpath = '/'.join([genomes_dir, genome_id, '%s.PATRIC.%s.tab' % (genome_id, type) ])
			if not os.path.isfile(inpath):
				print("%s not found" % os.path.basename(inpath) )
				continue

			infile = open(inpath)
			next(infile)
			for line in infile:
			
				# parse record
				r = line.rstrip().split('\t')
				genome_id = r[0]
				accession = r[2]
				seed_id = r[5]
				coord1 = r[9]
				coord2 = r[10]
				
				# skip genes missing identifier
				if seed_id == '': continue

				# format record
				gene_id = seed_id.split('|')[1] # remove fig| from identifier
				sid_long = 'accn|%s' % accession
				start = min(int(coord1), int(coord2))-1 # 0-based start coord
				end = max(int(coord1), int(coord2))     # 1-based stop coord
				
				# write bedfile
				if gene_id in gene_to_pangene: # only record pangenes
					pangene_id = gene_to_pangene[gene_id]
					outfile.write('\t'.join([sid_long, str(start), str(end), gene_id, pangene_id])+'\n')
				else: # some genes not included in pangenome; not sure why this would happen
					print("Gene %s in genome %s missing from pangenome %s" % (gene_id, genome_id, cluster_id) )





