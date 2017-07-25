# Tutorial

This step-by-step tutorial will walk you through downloading, installing, and running MIDAS

## Download and installation
MIDAS is written in Python and runs on Linux and OSX  

Download the software:  
`git clone https://github.com/snayfach/MIDAS`  
[read more...](install.md)

Install python dependencies as needed:  
`python MIDAS/setup.py install`  
[read more...](install.md)

Download & unpack reference database:  
[http://lighthouse.ucsf.edu/MIDAS/midas_db_v1.2.tar.gz](http://lighthouse.ucsf.edu/MIDAS/midas_db_v1.2.tar.gz)  
`tar -zxvf midas_db_v1.2.tar.gz`  
[read more...](ref_db.md)

Update your environment:  
`export PYTHONPATH=$PYTHONPATH:MIDAS`  
`export PATH=$PATH:MIDAS/scripts`   
`export MIDAS_DB=midas_db_v1.2`  

Optionally, download & unpack example dataset:  
[http://lighthouse.ucsf.edu/MIDAS/example.tar.gz](http://lighthouse.ucsf.edu/MIDAS/example.tar.gz)  
`tar -zxvf example.tar.gz`


## Run MIDAS

Running MIDAS can be conceptually broken down in three steps:  
1) run MIDAS per sample: `run_midas.py [species, genes, snps]`  
* run_midas.py species should be run prior to genes or snps  
* run_midas.py genes and snps are independent of eachother and do not need to both be run  

2) merge results across samples: `merge_midas.py [species, genes, snps]`  
* creates output files that can facilitate comparative analysis across samples and species/genes/snps  

3) analyze results: `gene_diversity.py, snp_sharing.py, etc.`  
* these are the scripts currently in `/path/to/MIDAS/scripts/analysis`  
* at the moment they experimental but will be available directly on `/path/to/MIDAS/scripts` in the next update

### 1. Run MIDAS per-sample

First, create a new directory to store per-sample outputs:  
`mkdir midas_output`  

The basic command to run midas per sample is:  
 `run_midas.py {species, genes, snps} outdir [options]`
 
You can use the `-h` flag to get more info on any of the commands:  
`run_midas.py -h`  
`run_midas.py species -h`  
`run_midas.py genes -h`  
`run_midas.py snps -h`   

#### A) Profile species abundances
`run_midas.py species midas_output/sample_1 -1 example/sample_1.fq.gz`  
`run_midas.py species midas_output/sample_2 -1 example/sample_2.fq.gz`

* Output files contain estimated species abundances for individual samples
* This enables automatically profiling strain-level variation of all species in downstream modules 
* This step can be skipped by directly specifying species ids (`--species_id`) when running `run_midas.py genes or snps`
* Description of output files can be found in `midas_output/sample_1/species/readme.txt`
* [read more...](species.md)

#### B) Profile strain-level gene content of abundant species  
`run_midas.py genes midas_output/sample_1 -1 example/sample_1.fq.gz`  
`run_midas.py genes midas_output/sample_2 -1 example/sample_2.fq.gz` 

* Output files contain predicted gene +/- and estimated gene copy # for all abundant species in individual samples
* Requires output from (1-A) unless you directly specify species using `--species_id` 
* Description of output files can be found in `midas_output/sample_1/genes/readme.txt`
* [read more...](cnvs.md)

#### C) Profile strain-level nucleotide variants of abundant species
`run_midas.py snps midas_output/sample_1 -1 example/sample_1.fq.gz`  
`run_midas.py snps midas_output/sample_2 -1 example/sample_2.fq.gz` 

* Output files contain genome-wide nucleotide variation statistics for all abundant species in individual samples
* Requires output from (1-A) unless you directly specify species using `--species_id`  
* Description of output files can be found in `midas_output/sample_1/snps/readme.txt`
* [read more...](snvs.md)


### 2. Merge MIDAS results across samples

The basic command usage of merge scripts is:  
 `merge_midas.py {species, genes, snps} outdir -i input -t intype [options]`

* `-i` indicates the sample directories output by run_midas.py  
* `-t` indicates the input type
* The output of merge commands are matrix files which facilitate comparative, quantitative analyses
* The merge command is necessary in order to identify SNPs

You can use the `-h` flag to get more info on any of the commands:  
`merge_midas.py -h`  
`merge_midas.py species -h`  
`merge_midas.py genes -h`  
`merge_midas.py snps -h`   


#### A) Merge species abundance across samples 

`merge_midas.py species merged_species -i midas_output/sample_1,midas_output/sample_2 -t list`     

* Requires that you've already run (1-A)  
* Description of output files can be found in `merged_species/readme.txt`
* [read more...](merge_species.md)

#### B) Merge strain-level pan-genome results across samples 
`merge_midas.py genes -i midas_output/sample_1,midas_output/sample_2 -t list`  

* Requires that you've already run (1-B)  
* Description of output files can be found in `merged_genes/Bacteroides_vulgatus_57955/readme.txt`
* [read more...](merge_cnvs.md)

#### C) Merge strain-level nucleotide variant results across samples  
`merge_midas.py snps ./merged_snps -i midas_output/sample_1,midas_output/sample_2 -t list`    

* This command performs core-genome, pooled SNP calling
* Requires that you've already run (1-C)  
* Description of output files can be found in `merged_snps/Bacteroides_vulgatus_57955/readme.txt`
* [read more...](merge_snps.md)
