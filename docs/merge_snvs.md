## Overview
Description: merge single-nucleotide variants for an individual species across samples  
Input: list of sample directories  
Output: core-genome SNPs, SNP annotations, SNP allele frequency matrix, SNP alternate alleles, SNP depth, core-genome consensus sequences, and a phylogenetic tree  


## Usage
```
optional arguments:
  -h, --help            show this help message and exit
  --threads INT         number of CPUs to use for merging files (1)
                        increases speed when merging across many samples

Input/Output:
  -i INPUT              input to sample directories output by run_midas.py genes
                        see '-t' for details
  -t {list,file,dir}    'list': -i is a comma-separated list of paths to sample directories (ex: /sample1,/sample2)
                        'dir': -i is a  directory containing all samples (ex: /samples_dir)
                        'file': -i is a file containing paths to sample directories (ex: sample_paths.txt)
  -o OUTDIR             output directory

Species filters (select subset of species from INPUT):
  --min_samples INT     all species with >= MIN_SAMPLES (1)
  --species_id CHAR     comma-separated list of species ids
                        a list of prevalent species can be obtained by running 'merge_midas.py species'
                        a map of species ids to species names can be found in 'ref_db/annotations.txt'
  --max_species INT     maximum number of species to merge (use all)

Sample filters (select subset of samples from INPUT):
  --sample_depth FLOAT  minimum average read depth per sample (5.0)
  --fract_cov FLOAT     fraction of reference sites covered by at least 1 read (0.4)
  --max_samples INT     maximum number of samples to process.
                        useful for quick tests (use all)

Site filters (select subset of genomic sites from INPUT):
  --site_depth INT      minimum number of mapped reads per site.
                        a high value like 20 will result in accurate allele frequencies, but may discard many sites.
                        a low value like 1 will retain many sites but may not result in accurate allele frequencies (3)
  --site_prev FLOAT     site has at least <site_depth> coverage in at least <site_prev> proportion of samples.
                        a value of 1.0 will select sites that have sufficent coverage in all samples.
                        a value of 0.0 will select all sites, including those with low coverage in many samples 
                        NAs recorded for included sites with less than <site_depth> in a sample (0.95)
  --site_maf FLOAT      minimum minor allele frequency of site across samples.
                        setting this to zero (default) will keep invariant sites across samples.
                        setting this above zero (e.g. 0.01, 0.02, 0.05) will only keep common variants
  --max_sites INT       maximum number of sites to include in output.
                        useful for quick tests (use all)
```

## Examples
1) Merge results for all species. Provide list of paths to sample directories:
`merge_midas.py snps -o outdir -i sample_1,sample_2 -t list`

2) Merge results for one species (id=57955):
`merge_midas.py snps --species_id 57955 -o outdir/57955 -i sample_1,sample_2 -t list`

3) Only use samples with >15x average depth and only use sites covered by >=10 reads in at least >=95% of samples:
`merge_midas.py snps -o outdir -i /path/to/samples -t dir --sample_depth 15 --site_depth 10 --site_prev 0.95`

4) Run a quick test:
`merge_midas.py snps -o outdir -i /path/to/samples -t dir --max_species 1 --max_samples 10 --max_sites 1000`

## Outputs
This module generates the following output files: 

* **{species_id}.snps.list**: lists of SNPs that passed QC
* **{species_id}.snps.info**: info on SNPs that passed QC
* **{species_id}.snps.ref_freq**: reference allele frequency matrix (snps x samples)
* **{species_id}.snps.depth**: read depth matrix (snps x samples)
* **{species_id}.snps.fasta**: fasta file containing consensus sequences for species in each sample 
* **{species_id}.tree**: maximum-likelihood phylogenetic tree built from consensus sequences

## Memory usage  
* Memory usage will depend on the number of sites and samples in the resulting output files
* Phylogenetic tree building is the most memory-intensive step. For more info, see: http://www.microbesonline.org/fasttree/#GenomeWide
