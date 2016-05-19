# Tutorial

This step-by-step tutorial will walk you through downloading, installing, and running MIDAS

## Download and installation
MIDAS is written in Python and runs on Linux and OSX  
[read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/requires.md)  

Download the software:  
`git clone https://github.com/snayfach/MIDAS`  
[read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/install.md)

Install python dependencies:  
`python MIDAS/setup.py install`  
[read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/install.md)

Update your PATH:  
`export PATH=$PATH:MIDAS/scripts` 

Download the reference database:  
`python MIDAS/scripts/download_ref_db.py`   
[read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/ref_db.md)  

## Run MIDAS

Running MIDAS can be conceptually broken down in three steps:  
1) run MIDAS per sample: `run_midas.py [species, genes, snps]`  
2) merge results across samples: `merge_species.py, merge_genes.py, merge_snps.py`  
3) analyze results: `genome_diversity.py, gene_diversity.py, snp_sharing.py, etc.`  

First, move to the example directory and create a new directory to store per-sample output for species, genes, and SNPs:  
`cd MIDAS/example`  
`mkdir samples`  

###Run MIDAS per-sample

Basic command usage (see below for details):  
 `run_midas.py [species, genes, snps] outdir [options]`

Note that `outdir` should be named with a unique sample identifier

For options, examples, and more info, use the `-h` flag:  
`run_midas.py -h`  
`run_midas.py species -h`  
`run_midas.py genes -h`  
`run_midas.py snps -h`   

Each script should be run on each metagenomic sample: 

**Profile species abundances**  
`run_midas.py species samples/sample_1 -1 sample_1.fq.gz`  
`run_midas.py species samples/sample_2 -1 sample_2.fq.gz`  

* This enables automatically profiling strain-level variation of all species in downstream modules. 
* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/species.md)

**Profile strain-level gene content of abundant species**  
`run_midas.py genes samples/sample_1 -1 sample_1.fq.gz`  
`run_midas.py genes samples/sample_2 -1 sample_2.fq.gz` 

* Output files contain estimated pan-genome gene copy numbers for all abundant species 
* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/cnvs.md)

**Profile strain-level nucleotide variants of abundant species**  
`run_midas.py snps samples/sample_1 -1 sample_1.fq.gz`  
`run_midas.py snps samples/sample_2 -1 sample_2.fq.gz` 

* Output files contain genome-wide nucleotide variation statistics for all abundant species 
* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/snvs.md)


###Merge MIDAS results across samples

Basic command usage (see below for details):  
 `merge_midas.py [species, genes, snps] -i input -t intype -o example`

* `-i` indicates the sample directories output by run_midas.py  
* `-t` indicates the input type 
* `-o` is the basename for output files 

For options, examples, and more info, use the `-h` flag:  
`merge_midas.py -h`  
`merge_midas.py species -h`  
`merge_midas.py genes -h`  
`merge_midas.py snps -h`   

**Merge species abundance across samples**  
`merge_species.py -i samples/sample_1,samples/sample_2 -t list -o example`  

* `-i` is a list of sample directories  and `-t list` indicates that `-i` is a list of paths  
* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/merge_species.md)

**Merge strain-level pan-genome results across samples**  
`merge_genes.py -s 57955 -i samples/sample_1,samples/sample_2 -t list -o genes/57955`

* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/merge_cnvs.md)

**Merge strain-level nucleotide variant results across samples**  
`merge_snps.py -s 57955 -i samples/sample_1,samples/sample_2 -t list -o snps/57955`

* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/merge_snvs.md)
