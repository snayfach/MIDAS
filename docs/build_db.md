# Build custom MIDAS database

This script will allow you to build a MIDAS database using your own genomes or metagenomic assemblies for one or more species. 

## How it works

This script uses user-defined species and user-supplied genomes to build database files for running MIDAS. This will allow you to estimate the abundance of these species in a shotgun metagenome and quantify gene content and SNPs for species with sufficient data.

First, a pan-genome is built for each species. A pan-genome is defined as the set of non-redundant genes across all genomes for each species. The program USEARCH (http://www.drive5.com/usearch) is used for clustering genes. The default clustering identity is 95%

Next, a representative-genome database is built. Representative genomes are used for mapping reads and calling SNPs. The user needs to define which genome they want to use as the representative for each species. 

Finally, a marker genes database is built. Marker genes are defined as universal, single-copy gene families. These are genes that occur once per genome and in all genomes (of bacteria). MIDAS uses a set of 15 of such gene families. These are a subset of the PhyEco gene families described here: http://dx.doi.org/10.1371/journal.pone.0077033. To identify these genes, HMMER (http://hmmer.org) is used to scan each species' pan-genome. Once identified, a HS-BLASTN (http://dx.doi.org/10.1093/nar/gkv784) database is built for mapping short reads.

## Requirements
As with all scripts, MIDAS and its dependencies will need to be installed.  
Additionally, you need to have the following command-line tools installed:

* USEARCH: http://www.drive5.com/usearch/download.html
* HMMER3: http://hmmer.org

After installing them, you need to add their installation directories to your PATH:  
`export PATH=$PATH:/usearch/installation/directory`  
`export PATH=$PATH:/hmmer/installation/directory`

You should be able to call these programs from your command line:  
`usearch`  
`hmmsearch`

## Usage

`build_midas_db.py indir mapfile outdir [options]`

### Input/output files:

<b>indir PATH</b>
Path to directory of genomes. Each subdirectory should be named with a genome identifier. Each subdirectory should contain the following files:
	
\<genome_id>.fna
Genomic DNA sequence in FASTA format    

\<genome_id>.ffn
Gene DNA sequences in FASTA format 
 
\<genome_id>.features
File specificy genomic coordinates of genes.
This file should be tab-delimited with a header and 5 fields:  

* scaffold_id (CHAR)
* gene_id (CHAR)
* start (INT)
* end (INT)
* strand ('+' or '-')	 
	
<b>mapfile PATH</b>  
Path to mapping file that specifies which genomes belonging to the same species.  
The file should be tab-delimited file with a header and 3 fields:  

* genome_id (CHAR): corresponds to subdirectory within INDIR 
* species_id (CHAR): : species identifier for genome_id
* rep_genome (0 or 1): indicator if genome_id should be used for SNP calling

<b>outdir PATH</b>  
Directory to store MIDAS database

### Options:

<b>--type {marker,gene,genome,all}</b>  
Type of database to build (default=all)  

* <i>marker</i>: build database of marker genes for metagenomic species profiling. enables you to use the command "run_midas.py species"  
* <i>gene</i>: build database of genes for metagenomic pan-genome profiling   enables you to use the command "run_midas.py genes"  
* <i>genome</i>: build database of representative genomes for metagenomic SNP profiling  enables you to use the command "run_midas.py snps"  
* <i>all</i>: build all three types of databases sequentially  

<b>--threads INT</b>  
Number of threads to use

<b>--pid FLOAT</b>  
Percent identity threshold for clustering genes (default=0.95)

<b>--iter_size INT</b>  
Maximum number of genomes to process at one time to prevent exceeding USEARCH's 4G memory limit (default=500).  
If the number of genomes per species exceeds this parameter, USEARCH will be run in batches of <iter_size>. After all batches have completed, the resulting centroids from each batch are combined and clustered one final time.

<b>--max_species INT</b>  
Maximum number of species to process from input (default=use all).  
Useful for quick tests

<b>--max_genomes INT</b>  
Maximum number of genomes to process per species (default=use all).  
Useful for quick tests

<b>--compress</b>              
Compress output files with gzip
