## Overview
Description: merge single-nucleotide variants for an individual species across samples  
Input: list of sample directories  
Output: core-genome SNPs, SNP annotations, SNP allele frequency matrix, SNP alternate alleles, SNP depth, core-genome consensus sequences, and a phylogenetic tree  


## Usage
```
optional arguments:
  -h, --help            show this help message and exit
  --threads THREADS     Number of threads to use

Input/Output:
  -i INPUT              input to sample directories output by run_phylo_cnv.py genes
                        see '-t' for details
  -t {list,file,dir}    'list': -i is a comma-separated list of paths to sample directories (ex: /sample1,/sample2)
                        'dir': -i is a  directory containing all samples (ex: /samples_dir)
                        'file': -i is a file containing paths to sample directories (ex: sample_paths.txt)
  -s SPECIES_ID         species identifier
                        a list of prevalent species can be obtained by running 'merge_species.py'
                        a map of species ids to species names can be found in 'ref_db/annotations.txt'
  -o OUTDIR             output directory

Pipeline options (choose one or more; default=all):
  --snps                identify and store list of hq snps
  --freq                build allele frequency & depth matrixes
  --cons                generate fasta file of consensus sequences
  --tree                build phylogenetic tree

Sample filters (select subset of samples from INPUT):
  --sample_depth SAMPLE_DEPTH
                        minimum average read depth per sample (5.0)
  --fract_cov FRACT_COV
                        fraction of reference sites covered by at least 1 read (0.4)
  --max_samples MAX_SAMPLES
                        maximum number of samples to process.
                        useful for quick tests (use all)

Site filters (select subset of genomic sites from INPUT):
  --site_depth SITE_DEPTH
                        minimum number of mapped reads per site.
                        a high value like 20 will result in accurate allele frequencies, but may discard many sites.
                        a low value like 1 will retain many sites but may not result in accurate allele frequencies (3)
  --site_prev SITE_PREV
                        site has at least <site_depth> coverage in at least <site_prev> proportion of samples.
                        a value of 1.0 will select sites that have sufficent coverage in all samples.
                        a value of 0.0 will select all sites, including those with low coverage in many samples 
                        NAs recorded for included sites with less than <site_depth> in a sample (0.95)
  --max_sites MAX_SITES
                        maximum number of sites to include in output.
                        useful for quick tests (use all)
```

## Examples
Examples:
1) Merge results for species 57955. Provide list of paths to sample directories:  
`merge_snps.py -s 57955 -o outdir/57955 -i sample_1,sample_2 -t list`

2) Only use samples with >15x average depth and only use sites covered by >=10 reads in at least >=95% of samples:  
`merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --sample_depth 15 --site_depth 10 --site_prev 95`

3) Just identify core-genome sites and build matrixes; do not build consensus seqs or phylogenetic trees:  
`merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --snps --freq`

4) Run a quick test:  
`merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --max_samples 10 --max_sites 1000`

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
