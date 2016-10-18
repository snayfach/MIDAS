# Tutorial

This step-by-step tutorial will walk you through downloading, installing, and running MIDAS

## Download and installation
MIDAS is written in Python and runs on Linux and OSX  

Download the software:  
`git clone https://github.com/snayfach/MIDAS`  
[read more...] (install.md)

Install python dependencies as needed:  
`python MIDAS/setup.py install`  
[read more...] (install.md)

Download & unpack reference database:  
[http://lighthouse.ucsf.edu/MIDAS/midas_db_v1.1.tar.gz](http://lighthouse.ucsf.edu/MIDAS/midas_db_v1.1.tar.gz)  
`tar -zxvf midas_db_v1.1.tar.gz midas_db_v1.1`  
[read more...] (ref_db.md)

Update your environment:  
`export PYTHONPATH=$PYTHONPATH:MIDAS`  
`export PATH=$PATH:MIDAS/scripts`   
`export MIDAS_DB=midas_db_v1.1`  

Optionally, download & unpack example dataset:  
[http://lighthouse.ucsf.edu/MIDAS/example.tar.gz](http://lighthouse.ucsf.edu/MIDAS/example.tar.gz)  
`tar -zxvf example.tar.gz example`


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

###Run MIDAS per-sample

First, create a new directory to store per-sample outputs:  
`mkdir midas_output`  

The basic command to run midas per sample is:  
 `run_midas.py {species, genes, snps} outdir [options]`
 
You can use the `-h` flag to get more info on any of the commands:  
`run_midas.py -h`  
`run_midas.py species -h`  
`run_midas.py genes -h`  
`run_midas.py snps -h`   

#### 1) Profile species abundances
`run_midas.py species midas_output/sample_1 -1 example/sample_1.fq.gz`  
`run_midas.py species midas_output/sample_2 -1 example/sample_2.fq.gz`

* This enables automatically profiling strain-level variation of all species in downstream modules. 
* [read more...] (species.md)

#### 2) Profile strain-level gene content of abundant species  
`run_midas.py genes midas_output/sample_1 -1 example/sample_1.fq.gz`  
`run_midas.py genes midas_output/sample_2 -1 example/sample_2.fq.gz` 

* Requires that you've already run (1)
* Output files contain estimated pan-genome gene copy numbers for all abundant species 
* [read more...] (cnvs.md)

#### 3) Profile strain-level nucleotide variants of abundant species
`run_midas.py snps midas_output/sample_1 -1 example/sample_1.fq.gz`  
`run_midas.py snps midas_output/sample_2 -1 example/sample_2.fq.gz` 

* Requires that you've already run (1)  
* Output files contain genome-wide nucleotide variation statistics for all abundant species 
* [read more...] (snvs.md)


###Merge MIDAS results across samples

The basic command usage of merge scripts is:  
 `merge_midas.py {species, genes, snps} outdir -i input -t intype [options]`

* `-i` indicates the sample directories output by run_midas.py  
* `-t` indicates the input type

You can use the `-h` flag to get more info on any of the commands:  
`merge_midas.py -h`  
`merge_midas.py species -h`  
`merge_midas.py genes -h`  
`merge_midas.py snps -h`   


#### 1) Merge species abundance across samples 

`merge_midas.py species merged_species -i midas_out/sample_1,midas_out/sample_2 -t list`     

The results will be in the folder named `merged_species`, and it will contain four files:  
```
count_reads.txt  
coverage.txt  
relative_abundance.txt  
species_prevalence.txt  
```
Those files can be looked at with your favorite spreadsheet or text editor for more careful analysis.

* [read more...] (merge_species.md)

#### 2) Merge strain-level pan-genome results across samples 
`merge_midas.py genes merged_genes -i midas_out/sample_1,midas_out/sample_2 -t list`  

The results will be in a folder named `merged_genes`, which will contain a sub-folder for each species that satisfies the coverage criterion for inclusion.

You can look into the folder with the output:
`ls merged_genes`
and it should contain a folder named `Bacteroides_vulgatus_57955`.
When looking at that folder, `ls merged_genes/Bacteroides_vulgatus_57955`, it should contain:  
```
genes_copynum.txt
genes_depth.txt
genes_info.txt
genes_presabs.txt
genes_summary.txt
```
Those files can be looked at with your favorite spreadsheet or text editor for more careful analysis.

* [read more...] (merge_cnvs.md)

#### 3) Merge strain-level nucleotide variant results across samples  
`merge_midas.py snps ./merged_snps -i midas_out/sample_1,midas_out/sample_2 -t list`    

The results will be in a folder named `./merged_snps`, which will contain a sub-folder for each species that satisfies the coverage criterion for inclusion.
For these example files, the screen should show:

You can look into the folder with the output:
`ls merged_snps`
and it should contain a folder named `Bacteroides_vulgatus_57955`.
When looking at that folder, `ls merged_snps/Bacteroides_vulgatus_57955` , it should contain:  
```
snps_alt_allele.txt  
snps_depth.txt  
snps_info.txt  
snps_ref_freq.txt  
snps_summary.txt
```
Those files can be looked at with your favorite spreadsheet or text editor for more careful analysis.

* [read more...] (merge_snvs.md)
