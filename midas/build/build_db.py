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
		self.dir='%s/%s' % (dir, id)
		self.files = {}
		self.init_files()
		self.is_rep = None

	def init_files(self):
		if not os.path.isdir(self.dir):
			sys.exit("\nError: genome directory '%s' does not exist" % (self.dir))
		for type in ['fna', 'ffn', 'faa', 'genes']:
			inpath = '%s/%s.%s' % (self.dir, self.id, type)
			if os.path.isfile(inpath):
				self.files[type] = inpath
			else:
				error = ""
				error += "\nError: could not locate input file '%s/%s.%s'\n" % (self.dir, self.id, type)
				error += "\nYour genome should contain the following files:\n"
				error += "  %s/%s.fna   (FASTA of genome sequence)\n" % (self.dir, self.id)
				error += "  %s/%s.ffn   (FASTA of gene sequences)\n" % (self.dir, self.id)
				error += "  %s/%s.faa   (FASTA of protein sequences)\n" % (self.dir, self.id)
				error += "  %s/%s.genes (Genomic coordinates of genes on genome)\n" % (self.dir, self.id)
				sys.exit(error)

class Gene:
	def __init__(self, id):
		""" Instantiate Gene """
		self.id = id
		self.genome_id = None
		self.seq = None
		self.length = None
		self.cluster_id = {}
		self.centroid_id = {}

class Pangenome:
	def __init__(self, sp, outdir, ext):
		""" Instantiate Pangenome """
		self.dir = '%s/pan_genomes/%s' % (outdir, sp.id)
		self.tmp = '%s/temp' % self.dir
		self.species = sp
		self.genomes = list(sp.genomes.values())
		self.stats = {}
		self.stats['genomes'] = len(self.genomes)
		self.count_genes = 0
		self.count_genes = 0
		try: os.makedirs(self.tmp)
		except: pass
	
	def store_genes(self, max_length):
		""" Store genes from all genomes """
		self.genes = {}
		self.stats['genes'] = 0
		for genome in self.genomes:
			for rec in Bio.SeqIO.parse(genome.files['ffn'], 'fasta'):
				if str(rec.seq) == '' or str(rec.id) in ['', '|']:
					continue
				elif len(str(rec.seq)) >= max_length:
					continue
				else:
					gene = Gene(rec.id)
					gene.genome_id = genome.id
					gene.seq = str(rec.seq).upper()
					gene.length = len(gene.seq)
					self.genes[gene.id] = gene
					self.stats['genes'] += 1

	def write_readme(self):
		""" Concatenate all genes from pangenome into sequence file """
		file = utility.iopen('%s/readme.txt' % self.dir, 'w')
		file.write("""
Description and statistics for pan-genome files

Summary Statistics
############

Genomes: %(genomes)s
Genes: %(genes)s
Gene clusters (99%% identity): %(centroids_99)s
Gene clusters (95%% identity): %(centroids_95)s
Gene clusters (90%% identity): %(centroids_90)s
Gene clusters (85%% identity): %(centroids_85)s
Gene clusters (80%% identity): %(centroids_80)s
Gene clusters (75%% identity): %(centroids_75)s
		
Output files
############
genes.ffn
  all genes from specified genomes
  
centroids.ffn
  gene sequences from 99%% identity gene clusters
  used for recruiting metagenomic reads
  
gene_info.txt
  information for all genes from genes.ffn
  the fields centroid_{99,95,90,95,80,75} indicate mappings between gene_id and gene clusters
""" % self.stats)
		file.close()

	def write_genes(self, resume):
		""" Concatenate all genes from pangenome into sequence file """
		ffn_path = '%s/genes.ffn' % self.dir
		if os.path.exists(ffn_path) and os.stat(ffn_path).st_size > 0 and resume:
			return
		file = utility.iopen('%s/genes.ffn' % self.dir, 'w')
		for gene in self.genes.values():
			file.write('>%s\n%s\n' % (gene.id, gene.seq))
		file.close()

	def cluster_genes(self, threads, resume):
		""" Cluster genes at 99% ID; Clustering centroids at lower %ID cutoffs """
		
		ffn_path = '%s/centroids.99.ffn' % self.tmp
		if not os.path.exists(ffn_path) or os.stat(ffn_path).st_size == 0 or not resume:
			self.uclust(
				genes='%s/genes.ffn' % self.dir,
				pid=0.99,
				centroids='%s/centroids.99.ffn' % self.tmp,
				clusters='%s/uclust.99.txt' % self.tmp,
				threads=threads)
		self.store_gene_info(pid=99)
		shutil.copy('%s/centroids.99.ffn' % self.tmp, '%s/centroids.ffn' % self.dir)
		
		for pid in [95, 90, 85, 80, 75]:
			ffn_path = '%s/centroids.%s.ffn' % (self.tmp, pid)
			if not os.path.exists(ffn_path) or os.stat(ffn_path).st_size == 0 or not resume:
				self.uclust(
					genes='%s/centroids.99.ffn' % self.tmp,
					pid=pid/100.0,
					centroids='%s/centroids.%s.ffn' % (self.tmp, pid),
					clusters='%s/uclust.%s.txt' % (self.tmp, pid),
					threads=threads)
			self.store_gene_info(pid)
		self.store_cluster_membership()

	def store_gene_info(self, pid):
		""" Parse UCLUST file and store mapping of gene_id to centroid_id at given %ID cutoff """
		self.stats['centroids_%s' % pid] = 0
		for r in self.parse_uclust('%s/uclust.%s.txt' % (self.tmp, pid)):
			if r['type'] == 'H':
				self.genes[r['gene_id']].cluster_id[pid] = r['cluster_id']
				self.genes[r['gene_id']].centroid_id[pid] = r['centroid_id']
			elif r['type'] == 'S':
				self.genes[r['gene_id']].cluster_id[pid] = r['cluster_id']
				self.genes[r['gene_id']].centroid_id[pid] = r['gene_id']
				self.stats['centroids_%s' % pid] += 1
			else:
				continue
			
	def store_cluster_membership(self):
		""" Map gene to 99% ID centroids at each clustering %ID cutoff """
		for gene in self.genes.values():
			gene.centroid_99 = gene.centroid_id[99]
			gene.centroid_95 = self.genes[gene.centroid_99].centroid_id[95]
			gene.centroid_90 = self.genes[gene.centroid_99].centroid_id[90]
			gene.centroid_85 = self.genes[gene.centroid_99].centroid_id[85]
			gene.centroid_80 = self.genes[gene.centroid_99].centroid_id[80]
			gene.centroid_75 = self.genes[gene.centroid_99].centroid_id[75]

	def write_gene_info(self):
		""" Record gene info in gene_info.txt """
		file = utility.iopen('%s/gene_info.txt' % self.dir, 'w')
		header = ['gene_id', 'genome_id', 'gene_length', 'centroid_99', 'centroid_95', 'centroid_90', 'centroid_85', 'centroid_80', 'centroid_75']
		file.write('\t'.join(header)+'\n')
		for gene_id in sorted(self.genes.keys()):
			g = self.genes[gene_id]
			values = [g.id, g.genome_id, g.length, g.centroid_99, g.centroid_95, g.centroid_90, g.centroid_85, g.centroid_80, g.centroid_75]
			file.write('\t'.join([str(_) for _ in values])+'\n')
		file.close()

	def clean_up(self):
		""" Remove temporary files """
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
		""" Run UCLUST from shell with specified arguments """
		command = "vsearch "
		command += "-cluster_fast %s " % genes
		command += "-id %s " % pid
		command += "-centroids %s " % centroids
		command += "-uc %s " % clusters
		command += "-threads %s " % threads
		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		utility.check_exit_code(process, command)

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
		if len(line.rstrip()) == 0: continue
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
				sp.genomes[genome_id].is_rep = True
			else:
				sp.genomes[genome_id].is_rep = False
		# store species
		if len(species) < args['max_species']:
			species[species_id] = sp
	# update # of genomes
	for sp in species.values():
		sp.ngenomes = len(sp.genomes)
		# make sure at least 1 rep genome/species
		if sp.rep_genome is None:
			sp.rep_genome = sp.genomes.keys()[0]
	return list(species.values())

def read_genomes(species):
	genomes = sum([list(sp.genomes.values()) for sp in species], [])
	return genomes

def build_repgenome_db(args, genomes, species):
	for sp in species:
		print("%s" % sp.id)
		outdir = '%s/rep_genomes/%s' % (args['outdir'], sp.id)
		if not os.path.isdir(outdir): os.makedirs(outdir)
		shutil.copy(sp.genomes[sp.rep_genome].files['genes'], '%s/genome.features' % outdir)
		#build_features_file(sp, fpath='%s/genome.features' % outdir)
		shutil.copy(sp.genomes[sp.rep_genome].files['fna'], '%s/genome.fna' % outdir)
		
def find_gene(gene, contigs):
	fwd_gene = str(gene).upper()
	rev_gene = str(gene.reverse_complement()).upper()	
	for id, contig in contigs:
		for seq, strand in [(fwd_gene, '+'), (rev_gene, '-')]:
			try: 
				start = contig.index(seq) + 1
				end = start + len(seq) - 1
				return (id, start, end, strand)
			except:
				continue
	sys.exit("Gene not found")
				
def build_features_file(sp, fpath):
	
	contigs = [] 
	with open(sp.genomes[sp.rep_genome].files['fna']) as f:
		for contig in Bio.SeqIO.parse(f, 'fasta'):
			contigs.append([contig.id, str(contig.seq).upper()])
	
	features = []
	with open(sp.genomes[sp.rep_genome].files['ffn']) as f:
		for gene in Bio.SeqIO.parse(f, 'fasta'):
			contig, start, end, strand = find_gene(gene.seq, contigs)
			features.append([gene.id, contig, start, end, strand])
			# prune list of contigs to increase speed
			# genes must be in sorted order
#			while True:
#				if contig != contigs[0][0]:
#					contigs = contigs[1:]
#				else:
#					break
					
	with open(fpath, 'w') as f:
		f.write('\t'.join(['gene_id', 'scaffold_id', 'start', 'end', 'strand'])+'\n')
		for r in features:
			f.write('\t'.join([str(_) for _ in r])+'\n')
		

def build_pangenome_db(args, species):
	for sp in species:
		print("%s" % sp.id)
		p = Pangenome(sp, outdir=args['outdir'], ext=args['compress'])
		if os.path.exists('%s/readme.txt' % p.dir) and args['resume']:
			print("  nothing to do")
			continue
		print("  catting genes")
		p.store_genes(args['max_length'])
		p.write_genes(args['resume'])
		print("  clustering genes")
		p.cluster_genes(args['threads'], args['resume'])
		print("  writing gene info")
		p.write_gene_info()
		print("  removing temporary files")
		p.clean_up()
		p.write_readme()

def write_species_info(args, species):
	outfile = utility.iopen('%s/species_info.txt' % args['outdir'], 'w')
	header = ['species_id', 'rep_genome', 'count_genomes']
	outfile.write('\t'.join(header)+'\n')
	for sp in species:
		values = [str(_) for _ in [sp.id, sp.rep_genome, sp.ngenomes]]
		outfile.write('\t'.join(values)+'\n')
		
def write_genome_info(args, species):
	outfile = utility.iopen('%s/genome_info.txt' % args['outdir'], 'w')
	header = ['genome_id', 'species_id', 'rep_genome']
	outfile.write('\t'.join(header)+'\n')
	for sp in species:
		for genome_id in sp.genomes:
			rep_genome = '1' if genome_id == sp.rep_genome else '0'
			values = [genome_id, sp.id, rep_genome]
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

def build_marker_db(args, genomes, species):
	marker_genes = MarkerGenes(args['outdir'])
	print("  searching marker gene HMMs vs pangenomes")
	for sp in species:
		for genome in sp.genomes.values():
			hmmpath = '%s/%s.hmmsearch' % (marker_genes.tmp, genome.id)
			if not os.path.exists(hmmpath) or os.stat(hmmpath).st_size == 0 or not args['resume']:
				marker_genes.hmmsearch(
					inpath=genome.files['faa'],
					outpath=hmmpath,
					threads=args['threads'])
			fna = marker_genes.parse_fasta(genome.files['ffn'])
			for h in marker_genes.find_hits(
					inpath='%s/%s.hmmsearch' % (marker_genes.tmp, genome.id),
					max_evalue=1e-5,
					min_cov=0.00):
				gene = fna[h['query']].upper()
				info = [sp.id, genome.id, h['query'], len(gene), h['target']]
				marker_genes.info.write('\t'.join([str(_) for _ in info])+'\n')
				if genome.is_rep:
					marker_genes.fasta.write('>'+h['query']+'\n'+gene+'\n')
	marker_genes.info.close()
	marker_genes.fasta.close()
	print("  building blast database")
	marker_genes.build_hsblastn_db(args['hs-blastn'], args['resume'])
	print("  writing mapping cutoffs file")
	marker_genes.build_mapping_cutoffs()
	print("  removing temporary files")
	shutil.rmtree(marker_genes.tmp)

class MarkerGenes:
	def __init__(self, dir):
		self.dir = '%s/marker_genes' % dir
		self.tmp = '%s/temp' % self.dir
		if not os.path.isdir(self.tmp): os.makedirs(self.tmp)
		self.fasta = open('%s/phyeco.fa' % self.dir, 'w')
		self.info = open('%s/phyeco.map' % self.dir, 'w')
		self.header = ['species_id', 'genome_id', 'gene_id', 'gene_length', 'marker_id']
		self.info.write('\t'.join(self.header)+'\n')

	def hmmsearch(self, inpath, outpath, threads):
		command = "hmmsearch --noali --cpu %s " % threads
		command += "--domtblout %s " % outpath
		command += "%s/%s " % (os.path.dirname(__file__), 'phyeco.hmm')
		command += "%s > /dev/null" % inpath
		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		utility.check_exit_code(process, command)

	def parse_hmmsearch(self, p_in):
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

	def find_hits(self, inpath, max_evalue, min_cov):
		hits = {}
		for r in self.parse_hmmsearch(inpath):
			if r['evalue'] > max_evalue:
				continue
			elif min(r['qcov'], r['tcov']) < min_cov:
				continue
			if r['target'] not in hits:
				hits[r['target']] = r
			elif r['evalue'] < hits[r['target']]['evalue']:
				hits[r['target']] = r
		return list(hits.values())

	def hsblastn_index(self, fasta):
		command = "hs-blastn index %s " % fasta
		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env={'PATH':sys.path})
		utility.check_exit_code(process, command)
	
	def parse_fasta(self, p_in):
		""" Return lookup of seq_id to sequence for PATRIC genes """
		seqs = {}
		infile = utility.iopen(p_in)
		for r in Bio.SeqIO.parse(infile, "fasta"):
			seqs[r.id] = str(r.seq).upper()
		infile.close()
		return seqs

	def build_hsblastn_db(self, hsblastn, resume):
		header = '%s/phyeco.fa.header' % self.dir
		if os.path.exists(header) and os.stat(header).st_size > 0 and resume:
			return
		command = "%s index " % hsblastn
		command += " %s/phyeco.fa " % self.dir
		process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		utility.check_exit_code(process, command)
	
	def build_mapping_cutoffs(self):
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
		outfile = open('%s/phyeco.mapping_cutoffs' % self.dir, 'w')
		for marker_id, cutoff in cutoffs.items():
			outfile.write(marker_id+'\t'+str(cutoff)+'\n')
		outfile.close()

def run_pipeline(args):
		
	print("Reading species & genome info")
	species = read_species(args)
	write_species_info(args, species)
	genomes = read_genomes(species)
	write_genome_info(args, species)
	
	print("\nBuilding pangenome database")
	build_pangenome_db(args, species)
				
	print("\nBuilding representative genome database")
	build_repgenome_db(args, genomes, species)

	print("\nBuilding marker genes database")
	build_marker_db(args, genomes, species)

	print("")
	if args['compress']:
		print("Compressing data\n")
		compress(args['outdir'])




