#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, subprocess, sys, shutil, gzip
from midas import utility
import Bio.SeqIO

class Species:
	def __init__(self, id):
		self.id=id
		self.genomes={}
		self.rep_genome=None

class Genome:
	def __init__(self, id, dir):
		self.id=id
		self.dir='%s/%s'%(dir,id)
		self.files = self.init_files()

	def init_files(self):
		files = {}
		if not os.path.isdir(self.dir):
			sys.exit("\nError: genome directory '%s' does not exist" % (self.dir))
		for type in ['fna', 'ffn', 'features']:
			for ext in ['', '.gz', '.bz2']:
				inpath = '%s/%s.%s%s' % (self.dir, self.id, type, ext)
				if os.path.isfile(inpath):
					files[type] = inpath
			if type not in files:
				error = ""
				error += "\nError: could not locate input file '%s/%s.%s(.gz|.bz2)'\n" % (self.dir, self.id, type)
				error += "\nYour genome should contain the following files:\n"
				error += "  %s/%s.fna (FASTA of genome sequence)\n" % (self.dir, self.id)
				error += "  %s/%s.ffn (FASTA of gene sequences)\n" % (self.dir, self.id)
				error += "  %s/%s.features (Genomic coordinates of genes)" % (self.dir, self.id)
				sys.exit(error)
		return files

def parse_mapping_file(args):
	infile = utility.iopen(args['mapfile'])
	fields = next(infile).rstrip('\n').split('\t')
	for field in ['genome_id', 'species_id']:
		if field not in fields:
			sys.exit("Error: mapping file '%s' has no field labeled '%s'" % (args['mapfile'], field))
	for field in fields:
		if field not in ['genome_id', 'species_id', 'rep_genome']:
			sys.exit("Error: mapping file '%s' has unknown field labeled '%s'" % (args['mapfile'], field))
	for line in infile:
		values = line.rstrip('\n').split('\t')
		record = dict([(f,v) for f,v in zip(fields, values)])
		if len(values) < len(fields):
			sys.exit("Error: mapping file '%s' has different number of fields per row" % args['mapfile'])
		if 'rep_genome' in fields and record['rep_genome'] not in ['0', '1']:
			sys.exit("Error: mapping file '%s' has unknown value '%s' for field 'rep_genome'" % (args['mapfile'], record['rep_genome']))
		yield record

def read_species(args):
	species = {}
	for record in parse_mapping_file(args):
		species_id = record['species_id']
		genome_id = record['genome_id']
		# fetch species
		sp = species[species_id] if species_id in species else Species(species_id)
		# update genomes
		if len(sp.genomes) < args['max_genomes']:
			sp.genomes[genome_id] = Genome(genome_id, args['indir'])
			if record['rep_genome'] == '1':
				sp.rep_genome = genome_id
		# store species
		if len(species) < args['max_species']:
			species[species_id] = sp
	# update # of genomes
	for sp in species.values():
		sp.ngenomes = len(sp.genomes)
		# make sure at least 1 rep genome/species
		if sp.rep_genome is None:
			sp.rep_genome = sp.genomes.keys()[0]
	return species.values()

def read_genomes(species):
	genomes = sum([sp.genomes.values() for sp in species], [])
	return genomes

def parse_hmmsearch(p_in):
	""" Parse HMMER domblout files. Return data-type formatted dictionary """
	f_in = utility.iopen(p_in)
	for line in f_in:
		if line[0] == '#': continue
		x = line.rstrip().split()
		query = x[0]
		target = x[3]
		evalue = float(x[12])
		qcov = (int(x[20]) - int(x[19]) + 1)/float(x[2])
		tcov = (int(x[16]) - int(x[15]) + 1)/float(x[5])
		yield {'query':query, 'target':target, 'evalue':evalue, 'qcov':qcov, 'tcov':tcov, 'qlen':int(x[2]), 'tlen':int(x[5])}

def find_hits(args, species, max_evalue, min_cov):
	inpath = "%s/marker_genes/temp/%s.hmmsearch" % (args['outdir'], species.id)
	hits = {}
	for r in parse_hmmsearch(inpath):
		if r['evalue'] > max_evalue:
			continue
		elif min(r['qcov'], r['tcov']) < min_cov:
			continue
		if r['target'] not in hits:
			hits[r['target']] = r
		elif r['evalue'] < hits[r['target']]['evalue']:
			hits[r['target']] = r
	return hits.values()

def hmmsearch(args, species):
	command = "hmmsearch --noali --cpu %s " % args['threads']
	command += "--domtblout %s/marker_genes/temp/%s.hmmsearch " % (args['outdir'], species.id)
	command += "%s/%s " % (os.path.dirname(__file__), 'phyeco.hmm')
	command += "%s/pan_genomes/%s/centroids.faa > /dev/null" % (args['outdir'], species.id)
	return command

def hsblastn_index(fasta):
	command = "hs-blastn index %s " % fasta
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env={'PATH':sys.path})
	utility.check_exit_code(process, command)

def parse_fasta(p_in):
	""" Return lookup of seq_id to sequence for PATRIC genes """
	seqs = {}
	infile = utility.iopen(p_in)
	for r in Bio.SeqIO.parse(infile, "fasta"):
		seqs[r.id] = str(r.seq).upper()
	infile.close()
	return seqs

def build_fasta_db(args, species):
	seqfile = utility.iopen('%s/marker_genes/phyeco.fa' % args['outdir'], 'w')
	mapfile = utility.iopen('%s/marker_genes/phyeco.map' % args['outdir'], 'w')
	mapfile.write('\t'.join(['species_id', 'gene_id', 'gene_length', 'marker_id'])+'\n') # add back genome id
	for index, sp in enumerate(species):
		path = '%s/pan_genomes/%s/centroids.ffn' % (args['outdir'], sp.id)
		fna = parse_fasta(path)
		for h in find_hits(args, sp, max_evalue=1e-5, min_cov=0.70):
			gene = fna[h['query']].upper()
			mapfile.write('%s\t%s\t%s\t%s\n' % (sp.id, h['query'], len(gene), h['target']))
			seqfile.write('>'+h['query']+'\n'+gene+'\n')
	seqfile.close()
	mapfile.close()

def build_marker_db(args, genomes, species):
	outdir = '%s/marker_genes/temp' % args['outdir']
	if not os.path.isdir(outdir): os.makedirs(outdir)
	print("1. Searching marker gene HMMs vs pangenomes")
	for sp in species:
		print("   species: %s" % sp.id)
		command = hmmsearch(args, sp)
		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		utility.check_exit_code(process, command)
	print("2. Building FASTA of homologs")
	build_fasta_db(args, species)
	print("3. Building HS-BLASTN database of homologs")
	command = "%s index %s/marker_genes/phyeco.fa " % (args['hs-blastn'], args['outdir'])
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	utility.check_exit_code(process, command)
	print("4. Writing mapping cutoffs file")
	build_mapping_cutoffs(args)

def build_mapping_cutoffs(args):
	cutoffs = {
		'B000032':95.50,
		'B000039':94.75,
		'B000041':98.00,
		'B000062':97.25,
		'B000063':96.00,
		'B000065':98.00,
		'B000071':95.25,
		'B000079':98.00,
		'B000080':95.25,
		'B000081':97.00,
		'B000082':95.25,
		'B000086':96.75,
		'B000096':96.75,
		'B000103':95.25,
		'B000114':94.50,
	}
	outfile = open('%s/marker_genes/phyeco.mapping_cutoffs' % args['outdir'], 'w')
	for marker_id, cutoff in cutoffs.items():
		outfile.write(marker_id+'\t'+str(cutoff)+'\n')
	outfile.close()

def build_repgenome_db(args, genomes, species):
	for sp in species:
		print("species: %s" % sp.id)
		outdir = '%s/rep_genomes/%s' % (args['outdir'], sp.id)
		if not os.path.isdir(outdir): os.makedirs(outdir)
		for ext in ['fna', 'features']:
			inpath = sp.genomes[sp.rep_genome].files[ext]
			outpath = '%s/genome.%s' % (outdir, ext)
			shutil.copy(inpath, outpath)

def build_pangenome_db(args, species):
	for sp in species:
		print("species: %s" % sp.id)
		pangenome = Pangenome(sp, args['outdir'], args['iter_size'], args['compress'])
		pangenome.cluster_genes(args['pid'], args['threads'])
		pangenome.translate()
		pangenome.record_info()
		pangenome.clean_up()

class Gene:
	def __init__(self, id):
		self.id = id

class GeneCluster:
	def __init__(self, id):
		self.id = id

class Pangenome:
	def __init__(self, sp, outdir, iter_size, ext):
		self.dir = '%s/pan_genomes/%s' % (outdir, sp.id)
		self.species = sp
		self.genomes = sp.genomes.values()
		self.ngenomes = len(self.genomes)
		self.iter_size = iter_size
		self.niters = 1 + (self.ngenomes - 1)/self.iter_size
			
	def batch_seqs(self, index, outdir, genomes):
		seqfile = utility.iopen('%s/genes.ffn' % outdir, 'w')
		mapfile = utility.iopen('%s/gene_to_genome.txt' % outdir, 'w')
		start = index * self.iter_size
		stop = start + self.iter_size
		for genome in self.genomes[start:stop]:
			for rec in Bio.SeqIO.parse(genome.files['ffn'], 'fasta'):
				if str(rec.seq) == '' or str(rec.id) in ['', '|']:
					continue
				else:
					seqfile.write('>%s\n%s\n' % (rec.id, str(rec.seq).upper()))
					mapfile.write('%s\t%s\n' % (rec.id, genome.id))
		seqfile.close()
		mapfile.close()

	def translate(self):
		infile = utility.iopen('%s/centroids.ffn' % self.dir)
		outfile = utility.iopen('%s/centroids.faa' % self.dir, 'w')
		for r in Bio.SeqIO.parse(infile, 'fasta'):
			if len(r.seq) % 3 == 0:
				outfile.write('>'+r.id+'\n'+str(r.seq.translate()).rstrip('*')+'\n')
		infile.close()
		outfile.close()

	def cluster_genes(self, pid, threads):
		print("  1. Clustering genes")
		# split genes into batches and cluster
		for index in range(self.niters):
			print("     iteration %s" % index)
			outdir = '%s/temp/iter_%s' % (self.dir, index) if self.niters > 1 else '%s/temp' % self.dir
			if not os.path.isdir(outdir): os.makedirs(outdir)
			self.batch_seqs(index, outdir, self.genomes)
			self.uclust('%s/genes.ffn' % outdir, pid, '%s/centroids.ffn' % outdir, '%s/uclust.txt' % outdir, threads)
		# if split into > 1 batch, combine batches and recluster
		if self.niters > 1:
			print("  2. Combining centroids across clustering iterations")
			seqfile = utility.iopen('%s/temp/genes.ffn' % self.dir, 'w')
			mapfile = utility.iopen('%s/temp/gene_to_genome.txt' % self.dir, 'w')
			for i in range(self.niters):
				indir = '%s/temp/iter_%s' % (self.dir, i)
				for line in utility.iopen('%s/centroids.ffn' % indir): seqfile.write(line)
				for line in utility.iopen('%s/gene_to_genome.txt' % indir): mapfile.write(line)
			seqfile.close()
			mapfile.close()
			print("  3. re-clustering combined centroids")
			self.uclust('%s/temp/genes.ffn' % self.dir, pid, '%s/temp/centroids.ffn' % self.dir, '%s/temp/uclust.txt' % self.dir, threads)
		# move centroids to final dest
		shutil.move('%s/temp/centroids.ffn' % self.dir, '%s/centroids.ffn' % self.dir)

	def map_gene_to_genome(self):
		genes = {}
		for line in utility.iopen('%s/temp/gene_to_genome.txt' % self.dir):
			gene_id, genome_id = line.rstrip('\n').split('\t')
			genes[gene_id] = Gene(gene_id)
			genes[gene_id].genome_id = genome_id
		return genes

	def store_gene_info(self):
		genes = self.map_gene_to_genome()
		for r in self.parse_uclust('%s/temp/uclust.txt' % self.dir):
			if r['type'] == 'H':
				genes[r['gene_id']].cluster_id = r['cluster_id']
				genes[r['gene_id']].centroid = '0'
				genes[r['gene_id']].length = int(r['size'])
			elif r['type'] == 'S':
				genes[r['gene_id']].cluster_id = r['cluster_id']
				genes[r['gene_id']].centroid = '1'
				genes[r['gene_id']].length = int(r['size'])
		return genes

	def store_cluster_info(self):
		clusters = {}
		for r in self.parse_uclust('%s/temp/uclust.txt' % self.dir):
			if r['type'] == 'C':
				clusters[r['cluster_id']] = GeneCluster(r['cluster_id'])
				clusters[r['cluster_id']].centroid_id = r['gene_id']
				clusters[r['cluster_id']].size = int(r['size'])
		return clusters

	def update_info(self):
		if self.niters == 1: return
		for index in range(self.niters):
			for r in self.parse_uclust('%s/temp/iter_%s/uclust.txt' % (self.dir, index)):
				if r['type'] == 'H':
					gene = self.genes[r['gene_id']]
					cluster_id = self.genes[r['centroid_id']].cluster_id
					gene.cluster_id = cluster_id
					gene.centroid = '0'
					gene.length = r['size']
					self.clusters[cluster_id].size += 1

	def write_cluster_info(self):
		outfile = utility.iopen('%s/cluster_info.txt' % self.dir, 'w')
		outfile.write('\t'.join(['cluster_id', 'size', 'centroid'])+'\n')
		for cluster_id in sorted(self.clusters.keys()):
			cluster = self.clusters[cluster_id]
			outfile.write('\t'.join([cluster.id, str(cluster.size), cluster.centroid_id])+'\n')
		outfile.close()

	def write_gene_info(self):
		outfile = utility.iopen('%s/gene_info.txt' % self.dir, 'w')
		outfile.write('\t'.join(['gene_id', 'genome_id', 'cluster_id', 'centroid', 'length'])+'\n')
		for gene_id in sorted(self.genes.keys()):
			gene = self.genes[gene_id]
			outfile.write('\t'.join([gene.id, gene.genome_id, gene.cluster_id, gene.centroid, str(gene.length)])+'\n')

	def record_info(self):
		self.genes = self.store_gene_info()
		self.clusters = self.store_cluster_info()
		self.update_info()
		self.write_cluster_info()
		self.write_gene_info()

	def clean_up(self):
		shutil.rmtree('%s/temp' % self.dir)

	def parse_uclust(self, inpath):
		""" Yield formatted records from UCLUST output file """
		# centroids are type == 'S'
		# non-centroids are type == 'H'
		# clusters are type == 'C'
		fields = ['type', 'cluster_id', 'size', 'pid', 'strand', 'skip1', 'skip2', 'skip3', 'gene_id', 'centroid_id']
		with utility.iopen(inpath) as infile:
			for index, line in enumerate(infile):
				values = line.rstrip('\n').split('\t')
				record = dict([(f,v) for f,v in zip(fields, values)])
				yield record

	def uclust(self, genes, pid, centroids, clusters, threads):
		command = "usearch "
		command += "-cluster_fast %s " % genes
		command += "-id %s " % pid
		command += "-centroids %s " % centroids
		command += "-uc %s " % clusters
		command += "-threads %s " % threads
		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		utility.check_exit_code(process, command)

def write_species_info(args, species):
	outfile = utility.iopen('%s/species_info.txt' % args['outdir'], 'w')
	header = ['species_id', 'rep_genome', 'count_genomes']
	outfile.write('\t'.join(header)+'\n')
	for sp in species:
		values = [str(_) for _ in [sp.id, sp.rep_genome, sp.ngenomes]]
		outfile.write('\t'.join(values)+'\n')

def compress(outdir):
	for module in ['pan_genomes', 'rep_genomes']:
		for species in os.listdir('%s/%s' % (outdir, module)):
			indir = '%s/%s/%s' % (outdir, module, species)
			for file in os.listdir(indir):
				inpath = '%s/%s' % (indir, file)
				if inpath.split('.')[-1] != 'gz':
					outfile = utility.iopen('%s/%s.gz' % (indir, file), 'w')
					for line in utility.iopen(inpath):
						outfile.write(line)
					outfile.close()
					os.remove(inpath)

def run_pipeline(args):
		
	species = read_species(args)
	genomes = read_genomes(species)
	
	write_species_info(args, species)
	
	print("Building pangenome database")
	print("=====================")
	build_pangenome_db(args, species)
				
	print("\nBuilding representative genome database")
	print("=====================")
	build_repgenome_db(args, genomes, species)

	print("\nBuilding marker genes database")
	print("=====================")
	build_marker_db(args, genomes, species)

	print("")
	if args['compress']:
		print("Compressing data\n")
		compress(args['outdir'])




