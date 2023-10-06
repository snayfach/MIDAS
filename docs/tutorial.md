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
[https://midasdb.pollard.gladstone.org/uhgg/midas\_db_v1.2.tar.gz](https://midasdb.pollard.gladstone.org/uhgg/midas_db_v1.2.tar.gz)  
`tar -zxvf midas_db_v1.2.tar.gz`  
[read more...](ref_db.md)

Update your environment:  
`export PYTHONPATH=$PYTHONPATH:MIDAS`  
`export PATH=$PATH:MIDAS/scripts`   
`export MIDAS_DB=midas_db_v1.2`  

Download & unpack example dataset:  
[https://midasdb.pollard.gladstone.org/uhgg/midas\_db_v1.2.tar.gz](https://midasdb.pollard.gladstone.org/uhgg/midas_db_v1.2.tar.gz)  
`tar -zxvf example.tar.gz`


## Run MIDAS

Running MIDAS can be conceptually broken down in three steps:  

<b> 1. Run MIDAS per sample:</b>  `run_midas.py [species, genes, snps]`  

* `run_midas.py species` should be run prior to `genes` or `snps`  
* If you know what species you're interested in, you can skip `run_midas.py species` and specify them using the `--species_id` flag with `run_midas.py genes` and `snps`
* `run_midas.py genes` and `snps` are independent of eachother and do not need to both be run  

<b> 2. Merge results across samples:</b> `merge_midas.py [species, genes, snps]`  

* Creates output matrices that can facilitate comparative analysis across samples and species/genes/snps  

<b> 3. Analyze results:</b> 

* A few scripts are provided for downstream analysis, which are listed on the [Table of Contents](https://github.com/snayfach/MIDAS/blob/dev/README.md)
* Otherwise, the output files can be analyzed and visualized using your favorite tools

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

#### A) [Profile species abundances](species.md)
`run_midas.py species midas_output/sample_1 -1 example/sample_1.fq.gz`  
`run_midas.py species midas_output/sample_2 -1 example/sample_2.fq.gz`

* Output files contain estimated species abundances for individual samples
* This enables automatically profiling strain-level variation of all species in downstream modules 
* This step can be skipped by directly specifying species ids (`--species_id`) when running `run_midas.py genes` or `snps`
* Description of output files can be found in `midas_output/<sample_id>/species/readme.txt`

#### B) [Profile strain-level gene content of abundant species](cnvs.md) 
`run_midas.py genes midas_output/sample_1 -1 example/sample_1.fq.gz`  
`run_midas.py genes midas_output/sample_2 -1 example/sample_2.fq.gz` 

* Output files contain predicted gene +/- and estimated gene copy # for all abundant species in individual samples
* Requires output from (1-A) unless you directly specify species using `--species_id` 
* Description of output files can be found in `midas_output/<sample_id>/genes/readme.txt`

#### C) [Profile strain-level nucleotide variants of abundant species](snvs.md)
`run_midas.py snps midas_output/sample_1 -1 example/sample_1.fq.gz`  
`run_midas.py snps midas_output/sample_2 -1 example/sample_2.fq.gz` 

* Output files contain genome-wide nucleotide variation statistics for all abundant species in individual samples
* Requires output from (1-A) unless you directly specify species using `--species_id`  
* Description of output files can be found in `midas_output/<sample_id>/snps/readme.txt`


### 2. Merge MIDAS results across samples

The basic command usage of merge scripts is:  
 `merge_midas.py {species, genes, snps} outdir -i input -t intype [options]`

* `-i` indicates the sample directories output by `run_midas.py`
* `-t` indicates the input type
* The output of merge commands are matrix files which facilitate comparative, quantitative analyses
* The merge command is necessary in order to perform pooled-sample, core-genome SNP calling

You can use the `-h` flag to get more info on any of the commands:  
`merge_midas.py -h`  
`merge_midas.py species -h`  
`merge_midas.py genes -h`  
`merge_midas.py snps -h`   


#### A) [Merge species abundance across samples](merge_species.md)

`merge_midas.py species merged_species -i midas_output/sample_1,midas_output/sample_2 -t list`     

* Requires that you've already run (1-A)  
* Description of output files can be found in `merged_species/readme.txt`

#### B) [Merge strain-level pan-genome results across samples](merge_cnvs.md)
`merge_midas.py genes -i midas_output/sample_1,midas_output/sample_2 -t list`  

* Requires that you've already run (1-B)  
* Description of output files can be found in `merged_genes/<species_id>/readme.txt`

#### C) [Merge strain-level nucleotide variant results across samples](merge_snps.md)  
`merge_midas.py snps ./merged_snps -i midas_output/sample_1,midas_output/sample_2 -t list`    

* This command performs pooled-sample, core-genome SNP calling
* Requires that you've already run (1-C)  
* Description of output files can be found in `merged_snps/<species_id>/readme.txt`
