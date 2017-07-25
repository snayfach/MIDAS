# Build custom MIDAS database

This script will allow you to build a MIDAS database using your own genomes or metagenomic assemblies for one or more species. 

## How it works

This script uses user-defined species and user-supplied genomes to build database files for running MIDAS. This will allow you to estimate the abundance of these species in a shotgun metagenome and quantify gene content and SNPs for species with sufficient data.

First, a pan-genome is built for each species. A pan-genome is defined as the set of non-redundant genes across all genomes for each species. The program [VSEARCH](https://github.com/torognes/vsearch) is first used to cluster genes at 99% percent identity and identify a centroid gene sequence from each cluster. These "centroids" are used for recruiting metagenomic reads. This is done to reduce the number of genes searched while maintaining maximum mapping sensitivity. Next, VSEARCH is used to cluster the centroids from 99% identity gene clusters further at 95, 90, 85, 80, and 75 percent identity. After mapping reads to 99% identity centroids with 'run_midas.py genes', reads can be optionally aggregated into gene clusters at any of the lower clustering thresholds with 'merge_midas.py genes'.

Next, a representative-genome database is built. Representative genomes are used for mapping reads and calling SNPs. The user simply needs to define which genome they want to use as the representative for each species.

Finally, a marker genes database is built. Marker genes are defined as universal, single-copy gene families. These are genes that occur once per genome and in all genomes (of bacteria). MIDAS uses a set of 15 of such gene families. These are a subset of the PhyEco gene families described here: http://dx.doi.org/10.1371/journal.pone.0077033. To identify these genes, HMMER (http://hmmer.org) is used to scan genes from the representative genome. Once identified, a HS-BLASTN (http://dx.doi.org/10.1093/nar/gkv784) database is built for mapping short reads.

## Requirements
As with all scripts, MIDAS and its dependencies will need to be installed.  
Additionally, you need to have the following command-line tools installed:

* VSEARCH: https://github.com/torognes/vsearch
* HMMER3: http://hmmer.org

After installing them, you need to add their installation directories to your PATH:  
`export PATH=$PATH:/vsearch/installation/directory`  
`export PATH=$PATH:/hmmer/installation/directory`

You should be able to call these programs from your command line:  
`vsearch`  
`hmmsearch`

## Usage

`build_midas_db.py indir mapfile outdir [options]`

### Input/output files:

<b>indir PATH</b>
Path to directory of genomes. Each subdirectory should be named with a genome identifier. Each subdirectory should contain the following files:
	
\<genome_id>.fna
Genomic DNA sequence in FASTA format    

\<genome_id>.faa
Protein sequences in FASTA format

\<genome_id>.ffn
Gene sequences in FASTA format

\<genome_id>.genes
Tab delimited file with genomic coordinates of genes. The file should be tab-delimited file with a header and the following fields:  

* gene_id (CHAR)
* scaffold_id (CHAR)
* start (INT)
* end (INT)
* strand (+ or -)
* gene_type (CDS or RNA)
 
<b>mapfile PATH</b>  
Path to mapping file that specifies which genomes belonging to the same species.  
The file should be tab-delimited file with a header and 3 fields:  

* genome_id (CHAR): corresponds to subdirectory within INDIR 
* species_id (CHAR): : species identifier for genome_id
* rep_genome (0 or 1): indicator if genome_id should be used for SNP calling

<b>outdir PATH</b>  
Directory to store MIDAS database

### Options:

<b>--threads INT</b>  
Number of threads to use

<b>--max_species INT</b>  
Maximum number of species to process from input (default=use all).  
Useful for quick tests

<b>--max_genomes INT</b>  
Maximum number of genomes to process per species (default=use all).  
Useful for quick tests

<b>--compress</b>              
Compress output files with gzip
