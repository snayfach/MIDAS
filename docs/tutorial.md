# Tutorial

This step-by-step tutorial will walk you through downloading, installing, and running MIDAS

## Download and installation
MIDAS is written in Python and runs on Linux and OSX  

Download the software:  
`git clone https://github.com/snayfach/MIDAS`  
[read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/install.md)

Install python dependencies as needed:  
`python MIDAS/setup.py install`  
[read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/install.md)

Update your environment:  
`export PYTHONPATH=$PYTHONPATH:MIDAS`  
`export PATH=$PATH:MIDAS/scripts` 

Download the reference database:  
`python MIDAS/scripts/download_ref_db.py`   
[read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/ref_db.md)  

## Run MIDAS

Running MIDAS can be conceptually broken down in three steps:  
1) run MIDAS per sample: `run_midas.py [species, genes, snps]`  
2) merge results across samples: `merge_midas.py [species, genes, snps]`  
3) analyze results: `genome_diversity.py, gene_diversity.py, snp_sharing.py, etc.`  

First, move to the example directory and create a new directory to store per-sample outputs:  
`cd MIDAS/example`  
`mkdir samples`  

###Run MIDAS per-sample

For options, examples, and more info, use the `-h` flag:  
`run_midas.py -h`  
`run_midas.py species -h`  
`run_midas.py genes -h`  
`run_midas.py snps -h`   

Basic command usage (see below for details):  
 `run_midas.py {species, genes, snps} outdir [options]`
 
**1) Profile species abundances**  
`run_midas.py species samples/sample_1 -1 sample_1.fq.gz`  
`run_midas.py species samples/sample_2 -1 sample_2.fq.gz`  

* This enables automatically profiling strain-level variation of all species in downstream modules. 
* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/species.md)

**2) Profile strain-level gene content of abundant species**  
`run_midas.py genes samples/sample_1 -1 sample_1.fq.gz`  
`run_midas.py genes samples/sample_2 -1 sample_2.fq.gz` 

* Requires that you've already run (1)
* Output files contain estimated pan-genome gene copy numbers for all abundant species 
* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/cnvs.md)

**3) Profile strain-level nucleotide variants of abundant species**  
`run_midas.py snps samples/sample_1 -1 sample_1.fq.gz`  
`run_midas.py snps samples/sample_2 -1 sample_2.fq.gz` 

* Requires that you've already run (1)  
* Output files contain genome-wide nucleotide variation statistics for all abundant species 
* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/snvs.md)


###Merge MIDAS results across samples

For options, examples, and more info, use the `-h` flag:  
`merge_midas.py -h`  
`merge_midas.py species -h`  
`merge_midas.py genes -h`  
`merge_midas.py snps -h`   

Basic command usage (see below for details):  
 `merge_midas.py {species, genes, snps} outdir -i input -t intype [options]`

* `-i` indicates the sample directories output by run_midas.py  
* `-t` indicates the input type

**1) Merge species abundance across samples**  
`merge_midas.py species ./merged_species -i samples/sample_1,samples/sample_2 -t list`

* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/merge_species.md)

**2) Merge strain-level pan-genome results across samples**  
`merge_midas.py genes ./merged_genes -i samples/sample_1,samples/sample_2 -t list`

* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/merge_cnvs.md)

**3) Merge strain-level nucleotide variant results across samples**  
`merge_midas.py snps ./merged_snps -i samples/sample_1,samples/sample_2 -t list`

* [read more...] (https://github.com/snayfach/MIDAS/blob/master/docs/merge_snvs.md)
