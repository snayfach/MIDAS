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

Update your environment:  
`export PYTHONPATH=$PYTHONPATH:MIDAS`  
`export PATH=$PATH:MIDAS/scripts` 

Download the reference database:  
`python MIDAS/scripts/download_ref_db.py`   
[read more...] (ref_db.md)

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
* [read more...] (species.md)

**2) Profile strain-level gene content of abundant species**  
`run_midas.py genes samples/sample_1 -1 sample_1.fq.gz`  
`run_midas.py genes samples/sample_2 -1 sample_2.fq.gz` 

* Requires that you've already run (1)
* Output files contain estimated pan-genome gene copy numbers for all abundant species 
* [read more...] (cnvs.md)

**3) Profile strain-level nucleotide variants of abundant species**  
`run_midas.py snps samples/sample_1 -1 sample_1.fq.gz`  
`run_midas.py snps samples/sample_2 -1 sample_2.fq.gz` 

* Requires that you've already run (1)  
* Output files contain genome-wide nucleotide variation statistics for all abundant species 
* [read more...] (snvs.md)


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
The results will be in the folder named `./merged_species`, and it will contain four files:
```
count_reads.txt
coverage.txt
relative_abundance.txt
species_prevalence.txt
```
Those files can be looked at with your favorite spreadsheet or text editor for more careful analysis.

* [read more...] (merge_species.md)

**2) Merge strain-level pan-genome results across samples**  
`merge_midas.py genes ./merged_genes -i samples/sample_1,samples/sample_2 -t list`  
The results will be in a folder named `./merged_genes`, which will contain a sub-folder for each species that satisfies the coverage criterion for inclusion.
For these example files, the screen should show:
```
MIDAS: Metagenomic Intra-species Diversity Analysis System
version 1.0.0; github.com/snayfach/MIDAS
Copyright (C) 2015-2016 Stephen Nayfach
Freely distributed under the GNU General Public License (GPLv3)

===========Parameters===========
Script: merge_midas.py genes
Input: samples/sample_1,samples/sample_2
Input type: list
Species selection criteria:
  >= 1 high-coverage samples per species
Sample selection criteria:
  >=1.0 average read depth across detected genes
Gene quantification criterea:
  present (1): genes with copy number >=0.35
  absent (0): genes with copy number <0.35
  cluster genes at 95 percent identity

Identifying species
  found 1 species with sufficient high-coverage samples

Merging: Bacteroides vulgatus (id:57955) for 2 samples
  building pangenome matrices
    22080 genes from 16 genomes
    clustered into 16268 families at 95 percent id
  writing gene info file
  writing summary statistics
```
You can look into the folder with the output:
`ls merged_genes`
and it should contain a folder named `57955`.
When looking at that folder, `ls merged_genes/57955` , it should contain:
```
genes_copynum.txt
genes_depth.txt
genes_info.txt
genes_presabs.txt
genes_summary.txt
```
Those files can be looked at with your favorite spreadsheet or text editor for more careful analysis.

* [read more...] (merge_cnvs.md)

**3) Merge strain-level nucleotide variant results across samples**  
`merge_midas.py snps ./merged_snps -i samples/sample_1,samples/sample_2 -t list`  
The results will be in a folder named `./merged_snps`, which will contain a sub-folder for each species that satisfies the coverage criterion for inclusion.
For these example files, the screen should show:
```
MIDAS: Metagenomic Intra-species Diversity Analysis System
version 1.0.0; github.com/snayfach/MIDAS
Copyright (C) 2015-2016 Stephen Nayfach
Freely distributed under the GNU General Public License (GPLv3)

===========Parameters===========
Script: merge_midas.py snps
Input: samples/sample_1,samples/sample_2
Input type: list
Output directory: ./merged_snps
Species identifier: None
Number of CPUs to use: 1
Sample selection criteria:
  keep samples with >=5.0 average coverage across reference genome
  keep samples where >=40.0 percent of reference genome has non-zero coverage
Site selection criteria:
  site must be covered by at least 3 reads across 95.0 percent of samples

Identifying species
  found 1 species with sufficient high-coverage samples

Merging: Bacteroides vulgatus (id:57955) for 1 samples
  merging per-sample statistics
  merging per-site statistics
  extracting and annotating specified sites
  removing temporary files
  ```

You can look into the folder with the output:
`ls merged_snps`
and it should contain a folder named `57955`.
When looking at that folder, `ls merged_snps/57955` , it should contain:  
```
snps_alt_allele.txt  
snps_depth.txt  
snps_info.txt  
snps_ref_freq.txt  
snps_summary.txt
```
Those files can be looked at with your favorite spreadsheet or text editor for more careful analysis.

* [read more...] (merge_snvs.md)
