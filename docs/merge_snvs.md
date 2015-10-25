## Overview
Merge SNPs across samples for a given species

## Usage
```
usage: merge_snps.py [options]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         verbose

Input/Output:
  -i IN, --indir IN     input directory
  -o OUT, --outdir OUT  output directory
  -g GENOME_CLUSTER, --genome_cluster GENOME_CLUSTER
                        genome cluster id
  -m MATRIX, --matrix MATRIX
                        reference SNP matrix

Pipeline options (choose one or more; default=all):
  --snps                identify and store list of hq snps
  --freq                build allele frequency & depth matrixes
  --cons                generate fasta file of consensus sequences
  --tree                build phylogenetic tree

Sample filters:
  --sample_list SAMPLE_LIST
                        file of sample ids to include; each line should
                        contain one id
  --sample_depth SAMPLE_DEPTH
                        min average read depth per sample (2.0)
  --ref_coverage REF_COVERAGE
                        min coverage of reference genome per sample (0.4)

SNP filters:
  --snp_prev MIN_PREV   min fraction of samples that contain SNP (1.0)
  --snp_depth SNP_DEPTH
                        min # of reads supporting SNP per sample (1)
  --no_fixed            exclude SNPs with the same consensus allele across
                        samples (False)
  --max_snps MAX_SNPS   only use <= MAX_SNPS (use all)
```

## Examples
Use default parameters and run entire pipeline:  
`merge_snps.py -i snvs -o B_vulgatus -g 57955`  

* `snvs` is a directory that contains results from `run_phylo_cnv snvs`
* the name of each subdirectory should correspond to a sample_id
* `57955` is the species identifier for Bacteroides vulgatus

Filter out low-coverage samples (<10x):  
`merge_snps.py -i snvs -o B_vulgatus -g 57955 --sample_depth 10.0`  

Exclude variants that have the same consensus allele across all samples:  
`merge_snps.py -i snvs -o B_vulgatus -g 57955 --no_fixed`  

Exclude variants that have the same consensus allele across all samples:  
`merge_snps.py -i snvs -o B_vulgatus -g 57955 --no_fixed`  

Keep variants with missing data in 5% of samples:  
`merge_snps.py -i snvs -o B_vulgatus -g 57955 --snp_prev 0.95`  

Identify snps and build allele frequency matrix only:  
`merge_snps.py -i snvs -o B_vulgatus -g 57955 ---snps --freq`  

Identify consensus sequences and build phylogenetic tree only:  
`merge_snps.py -i snvs -o B_vulgatus -g 57955 ---cons --tree` 
 

## Outputs
This module generates the following output files: 

* **{species_id}.hq_snps**: lists of SNPs that passed QC
* **{species_id}.ref_freq**: reference allele frequency matrix (snps x samples)
* **{species_id}.depth**: read depth matrix (snps x samples)
* **{species_id}.fasta**: fasta file containing consensus sequences for species in each sample 
* **{species_id}.tree**: maximum-likelihood phylogenetic tree built from consensus sequences

## Next steps
[Functionally annotate SNVs] (https://github.com/snayfach/PhyloCNV/blob/master/docs/annotate_snvs.md)
[Estimate phylogenetic distance between metagenomes] (https://github.com/snayfach/PhyloCNV/blob/master/docs/pairwise_distances.md)  
[Estimate population genetic parameters]
