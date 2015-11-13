#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.2'

import argparse, sys, os, gzip, subprocess, tempfile

def print_copyright():
	print ("")
	print ("PhyloCNV: species abundance and strain-level genomic variation from metagenomes")
	print ("version %s; github.com/snayfach/PhyloCNV" % __version__)
	print ("Copyright (C) 2015 Stephen Nayfach")
	print ("Freely distributed under the GNU General Public License (GPLv3)")
	print ("")

def parse_arguments():
	""" Parse command line arguments """
	
	parser = argparse.ArgumentParser(
		usage='%s [options]' % os.path.basename(__file__),
		description="""Merge single-nucleotide variants for an individual species across samples. Outputs include: a list of high-quality sites, an allele frequency matrix, consensus sequences for each sample, and a phylogenetic tree"""
		)
	parser.add_argument('-v', '--verbose', action='store_true', default=False, help='verbose')
	
	io = parser.add_argument_group('Input/Output')
#	io.add_argument('-i', dest='indir', type=str, required=True,
#		help="""input directory.
#			each subdirectory should correspond to a different sample_id""")
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""input to results from 'run_phylo_cnv.py snvs'.
			see <intype> for details""")
	io.add_argument('-t', choices=['dir', 'file', 'list'], dest='intype', required=True,
		help="""input type.
			'dir': directory containing phylo_cnv results. each subdirectory should correspond to a different sample_id. for example: <directory>/<sample_id>
			'file': file containing paths to phylo_cnv results.	each line in the file should contain the full path to the results for a sample_id.
			'list': comma-separated list of paths to phylo_cnv results.
			""")
	io.add_argument('-s', dest='species_id', type=str, required=True,
		help="""species identifier. 
			a list of prevalent species can be obtained by running 'scripts/merge_species.py'.
			A map of species ids to species names can be found in 'ref_db/annotations.txt'""")
	io.add_argument('-o', dest='outdir', type=str, required=True,
		help="""output directory.
			output files: <outdir>/<species_id>.hq_snps, <outdir>/<species_id>.ref_freq, <outdir>/<species_id>.depth, <outdir>/<species_id>.fasta, <outdir>/<species_id>.tree""")
	#io.add_argument('-m', '--matrix', type=str,  help='reference SNP matrix')
	
	pipe = parser.add_argument_group('Pipeline options (choose one or more; default=all)')
	pipe.add_argument('--snps', default=False, action='store_true', help='identify and store list of hq snps')
	pipe.add_argument('--freq', default=False, action='store_true', help='build allele frequency & depth matrixes')
	pipe.add_argument('--cons', default=False, action='store_true', help='generate fasta file of consensus sequences')
	pipe.add_argument('--tree', default=False, action='store_true', help='build phylogenetic tree')
	
	sample = parser.add_argument_group('Sample filters\n(determine which samples are included in output)')
	sample.add_argument('--sample_depth', dest='sample_depth', type=float, default=5.0,
		help='minimum average read depth per sample (5.0)')
	sample.add_argument('--fract_cov', dest='fract_cov', type=float, default=0.4,
		help='fraction of reference sites covered by at least 1 read (0.4)')
	sample.add_argument('--max_samples', type=int,
		help='maximum number of samples to process. useful for quick tests (use all)')
				
	snps = parser.add_argument_group('Site filters\n(determine which reference-genome positions are included in output)')
	snps.add_argument('--site_depth', dest='site_depth', type=int, default=3,
		help="""minimum number of mapped reads per site. a high value like 20 will result in accurate allele frequencies, but may discard many sites. a low value like 1 will retain many sites but may not result in accurate allele frequencies (3)""")
	snps.add_argument('--site_prev', dest='min_prev', type=float, default=0.95,
		help="""site has at least <site_depth> coverage in at least <site_prev> proportion of samples.
			a value of 1.0 will select sites that have sufficent coverage in all samples.
			a value of 0.0 will select all sites, including those with low coverage in many samples  (0.95)""")
	snps.add_argument('--no_fixed', action='store_true', default=False,
		help="""exclude sites with the same consensus allele across samples. 
			this can be useful to reduce the size of the output datasets while retaining most of the information (False)""")
	snps.add_argument('--max_sites', dest='max_snps', type=int, default=float('Inf'),
		help="""maximum number of sites to include in output. useful for quick tests (use all)""")
					
	args = vars(parser.parse_args())
	args['db'] = '%s/ref_db/genome_clusters' % os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	check_args(args)
	return args

def check_args(args):
	""" Check validity of command line arguments """
	# turn on entire pipeline
	if not any([args['snps'], args['freq'], args['cons'], args['tree']]):
		args['snps'] = True
		args['freq'] = True
		args['cons'] = True
		args['tree'] = True

def print_arguments(args):
	print ("-------------------------------------------------------")
	print ("Merge SNVs Parameters")
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Species identifier: %s" % args['species_id'])
	print ("Pipeline options:")
	if args['snps']: print ("  identify sites")
	if args['freq']: print ("  build allele frequency matrix")
	if args['cons']: print ("  call consensus sequences")
	if args['tree']: print ("  build phylogenetic tree")
	print ("Sample selection criteria:")
	if args['sample_depth']:
		print ("  keep samples with >=%s average coverage across reference genome" % args['sample_depth'])
	if args['fract_cov']:
		print ("  keep samples where >=%s percent of reference genome has non-zero coverage" % (100*args['fract_cov']))
	if args['max_samples']:
		print ("  analyze up to %s samples" % args['max_samples'])
	print ("Site selection criteria:")
	print ("  site must be covered by at least %s reads across %s percent of samples" % (args['site_depth'], 100*args['min_prev']))
	if args['no_fixed']:
		print ("  exclude 'fixed' sites with the same consensus allele in all samples")
	else:
		print ("  keep 'fixed' sites with the same consensus allele in all samples")
	print ("-------------------------------------------------------")
	print ("")

def parse_snps_summary(inpath):
	""" Read in summary snp statistics for genome-clusters """
	snps_summary = {}
	infile = open(inpath)
	fields = next(infile).rstrip().split()
	for line in open(inpath):
		values = line.rstrip().split()
		rec = dict([(i,j) for i,j in zip(fields, values)])
		snps_summary[rec['cluster_id']] = rec
	return snps_summary

def list_samples(input, intype):
	if intype == 'dir':
		if not os.path.isdir(input):
			sys.exit("\nSpecified input directory does not exist:\n%s" % input)
		else:
			return([os.path.join(input, _) for _ in os.listdir(input)])
	elif intype == 'file':
		if not os.path.isfile(input):
			sys.exit("\nSpecified input file does not exist:\n%s" % input)
		else:
			return([x.rstrip() for x in open(input).readlines()])
	elif intype == 'list':
		return(input.split(','))

def identify_samples(args):
	""" Identify samples with sufficient depth for target genome-cluster """
	samples_dirs = []
	for sample_dir in list_samples(args['input'], args['intype']):
		# read in summary stats for sample across genome-clusters
		inpath = '/'.join([sample_dir, 'snps_summary_stats.txt'])
		sample_id = os.path.basename(sample_dir)
		if not os.path.isfile(inpath):
			print("  warning: no data for sample_id %s" % sample_id)
			continue
		snps_summary = parse_snps_summary(inpath)
		# check whether sample passes QC
		if args['species_id'] not in snps_summary:
			continue
		elif float(snps_summary[args['species_id']]['average_depth']) < args['sample_depth']:
			continue
		elif float(snps_summary[args['species_id']]['fraction_covered']) < args['fract_cov']:
			continue
		# sample passes qc
		else:
			samples_dirs.append(sample_dir)
			# only keep max_samples if specified
			if args['max_samples'] and len(samples_dirs) >= args['max_samples']:
				break
	if len(samples_dirs) == 0:
		sys.exit("Error: species_id failed to pass quality-control in all samples")
	return samples_dirs

def parse_snps(inpath):
	""" Yields formatted records from SNPs output """
	infile = gzip.open(inpath)
	next(infile)
	fields = ['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
			  'count_alleles', 'count_ref', 'count_alt', 'depth', 'ref_freq']
	for line in infile:
		yield dict([(i,j) for i,j in zip(fields, line.rstrip().split())])

def open_infiles(args, sample_dirs, sample_ids):
	""" Open SNP files for genome-cluster across all samples """
	infiles = {}
	for sample_dir, sample_id in zip(sample_dirs, sample_ids):
		inpath = '%s/snps/%s.snps.gz' % (sample_dir, args['species_id'])
		infiles[sample_id] = parse_snps(inpath)
	return infiles

def fetch_snp(snpfiles, sample_ids):
	""" Fetch SNP data across samples """
	snp = []
	for sample_id in sample_ids:
		snp.append(next(snpfiles[sample_id]))
	return snp

def is_fixed(snp, pass_qc):
	""" Determine if SNP has the same consensus allele across all samples that passed QC """
	cons_alleles = [snp_i['cons_allele'] for is_hq, snp_i in zip(pass_qc, snp) if is_hq]
	if all([allele == cons_alleles[0] for allele in cons_alleles]):
		return True
	else:
		return False

def snp_qc(snp, args):
	""" Determine whether nor not to keep snp for each sample """
	i = []
	for s in snp:
		if float(s['depth']) < args['site_depth']:
			i.append(False)
		else:
			i.append(True)
	return i

def read_hq_snps(args):
	""" Read in list of HQ snps from file"""
	snps = []
	inpath = '%s/%s.hq_snps' % (args['outdir'], args['species_id'])
	infile = open(inpath)
	next(infile)
	for line in infile:
		snps.append(line.rstrip().split()[0:2])
	return snps

def count_alt_alleles(snp):
	""" Get alternate allele(s) for snp """
	alleles = []
	for s in snp:
		if s['alt_allele'] != 'NA':
			alleles.append(s['alt_allele'])
	count_alleles = dict( (allele, alleles.count(allele)) for allele in set(alleles) )
	return(sorted(count_alleles.items(), reverse=True))

def prevalence(list):
	""" Count fraction of True in list """
	return float(sum(list))/len(list)

def write_ref_freq(args, sample_dirs, sample_ids):
	""" Write reference allele frequencies """
	# open output files
	outfiles = {}
	outfiles['freq_matrix'] = open('%s/%s.ref_freq' % (args['outdir'], args['species_id']), 'w')
	outfiles['depth_matrix'] = open('%s/%s.depth' % (args['outdir'], args['species_id']), 'w')
	outfiles['freq_matrix'].write('\t'.join(['ref_id', 'ref_pos'] + sample_ids)+'\n')
	outfiles['depth_matrix'].write('\t'.join(['ref_id', 'ref_pos'] + sample_ids)+'\n')
	snpfiles = open_infiles(args, sample_dirs, sample_ids) # open all input files
	hq_snps = read_hq_snps(args) # read in list of hq snps
	n = len(hq_snps)
	i = 0 # index position in hq_snps
	dummyfile = gzip.open('%s/snps/%s.snps.gz' % (sample_dirs[0], args['species_id']))
	next(dummyfile)
	for line in dummyfile: # outer loop: stop when eof reached
		snp = fetch_snp(snpfiles, sample_ids) # inner loop: get data for snp across all samples
		id = [snp[0]['ref_id'], snp[0]['ref_pos']]
		if i == n: # no more snps
			break
		elif id != hq_snps[i]: # snp not in list
			continue
		else: # snp in list
			outfiles['freq_matrix'].write('\t'.join(id+[str(s['ref_freq']) for s in snp])+'\n')
			outfiles['depth_matrix'].write('\t'.join(id+[str(s['depth']) for s in snp])+'\n')
			i += 1

def id_snps(args, sample_dirs, sample_ids):
	""" Identify SNPs that pass QC """
	snpfiles = open_infiles(args, sample_dirs, sample_ids) # open all input files
	outfile = open('%s/%s.hq_snps' % (args['outdir'], args['species_id']), 'w')
	outfile.write('\t'.join(['ref_id', 'ref_pos', 'alt_allele', 'count_alt'])+'\n')
	hq_snps = []
	dummyfile = gzip.open('%s/snps/%s.snps.gz' % (sample_dirs[0], args['species_id']))
	next(dummyfile)
	for line in dummyfile: # outer loop: stop when eof reached
		snp = fetch_snp(snpfiles, sample_ids) # inner loop: get data for snp across all samples
		pass_qc = snp_qc(snp, args) # True/False depending if sample passed/failed QC for SNP
		if prevalence(pass_qc) < args['min_prev']: # skip snps that do not occur in enough samples
			continue
		elif args['no_fixed'] and is_fixed(snp, pass_qc): # skip snps that are are fixed across samples that passed QC
			continue
		else:
			hq_snps.append([snp[0]['ref_id'], snp[0]['ref_pos']]) # keep track of snps that passed qc
			alt_allele_count = count_alt_alleles(snp) # count times each alternate allele occured
			n_alleles = len(alt_allele_count)
			alt_allele = alt_allele_count[0][0] if n_alleles > 0 else 'NA'
			record = [snp[0]['ref_id'], snp[0]['ref_pos'], alt_allele, str(n_alleles)]
			outfile.write('\t'.join(record)+'\n')
		if len(hq_snps) == args['max_snps']: # stop after reaching max_snps (for testing only)
			break
	outfile.close()

def write_consensus(args, sample_dirs, sample_ids):
	""" Write consensus sequences from samples """
	outfile = open('%s/%s.fasta' % (args['outdir'], args['species_id']), 'w')
	hq_snps = read_hq_snps(args) # get snp list
	n = len(hq_snps)
	for sample_dir, sample_id in zip(sample_dirs, sample_ids):
		i = 0 # index position in hq_snps
		outfile.write('>%s\n' % sample_id) # sequence header
		inpath = '%s/snps/%s.snps.gz' % (sample_dir, args['species_id'])
		for snp in parse_snps(inpath):
			if i == n: # no more snps in list
				break
			elif [snp['ref_id'], snp['ref_pos']] != hq_snps[i]: # snp not in list
				continue
			else: # snp in list
				outfile.write('-' if snp['cons_allele'] == 'NA' else snp['cons_allele']) # write base
				i += 1
		outfile.write('\n')

#	# add consensus seqs for reference
#	i = 0 # index position in snp list
#	outfile.write('>representative_genome\n') # sequence header
#	inpath = '%s/%s/snps/%s.snps.gz' % (args['indir'], samples[0], args['species_id'])
#	for snp in parse_snps(inpath):
#		if i == n: # no more snps in list
#			break
#		elif [snp['ref_id'], snp['ref_pos']] != hq_snps[i]: # snp not in list
#			continue
#		else: # snp in list
#			outfile.write(snp['ref_allele']) # write base
#			i += 1
#	outfile.write('\n')

def add_ref_alleles(args):
	""" Add partial genotypes from reference genomes """
	meta_snps = read_hq_snps(args) # get snp list
	n_snps = len(meta_snps) # length of snp list
	ref_snps = gzip.open(args['matrix']) # open file of reference genotypes
	max_genomes = 300 # maximum genomes to use from matrix
	# open tempfiles for reference genotypes
	tempfiles = [] #
	for line in ref_snps:
		genomes = line.rstrip().split()[2:max_genomes+1]
		for genome in genomes:
			file = open(tempfile.mkstemp()[1], 'w')
			file.write('>%s\n' % genome)
			tempfiles.append(file)
		break
	# write reference genotypes to temp files
	i = 0 # index position in snp list
	for line in ref_snps: # get partial genotypes for reference genomes
		ref_snp = line.rstrip().split()
		if i == n_snps: # no more snps in meta_snp list
			break
		elif ref_snp[0:2] != meta_snps[i]: # ref_snp not in meta_snp list
			continue
		else: # ref_snp in meta_snp list
			genotypes = ref_snp[2:max_genomes+1]
			for genotype, outfile in zip(genotypes, tempfiles):
				outfile.write(genotype)
			i += 1 # increment index
	# append reference genotypes to fasta
	outpath = '%s/%s.fasta' % (args['outdir'], args['species_id'])
	outfile = open(outpath, 'a')
	for f in tempfiles:
		f.close()
		f = open(f.name)
		outfile.write(f.read()+'\n')
		os.remove(f.name)

def build_tree(args):
	"""	Use FastTree to build phylogenetic tree of consensus sequences """
	inpath = '%s/%s.fasta' % (args['outdir'], args['species_id'])
	outpath = '%s/%s.tree' % (args['outdir'], args['species_id'])
	p = subprocess.Popen('FastTree -nt -boot 100 < %s > %s' % (inpath, outpath), shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	out, err = p.communicate()

if __name__ == '__main__':

	args = parse_arguments()
	if args['verbose']: print_copyright()
	if args['verbose']: print_arguments(args)
	if not os.path.isdir(args['outdir']): os.mkdir(args['outdir'])

	if args['verbose']: print("Identifying samples with species")
	sample_dirs = identify_samples(args) # id samples with sufficient depth
	sample_ids = [os.path.basename(_) for _ in sample_dirs]
	if args['verbose']: print("  %s samples with species" % len(sample_ids))

	if args['snps']:
		if args['verbose']: print("Identifying and writing hq snps")
		id_snps(args, sample_dirs, sample_ids)
	
	if args['freq']:
		if args['verbose']: print("Writing allele frequencies & depths")
		write_ref_freq(args, sample_dirs, sample_ids)

	if args['cons']:
		if args['verbose']: print("Writing consensus sequences")
		write_consensus(args, sample_dirs, sample_ids)

#	if args['matrix']:
#		if args['verbose']: print("Adding references sequences")
#		add_ref_alleles(args)

	if args['tree']:
		if args['verbose']: print("Building phylogenetic tree")
		build_tree(args)

			
