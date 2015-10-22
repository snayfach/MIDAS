## Overview
Merge SNPs across samples for a given species

## Usage
```
usage: merge_snps.py [options]

optional arguments:
  -h, --help            show this help message and exit

Input/Output:
  -i IN, --indir IN     input directory
  -o OUT, --outdir OUT  output directory
  -g GENOME_CLUSTER, --genome_cluster GENOME_CLUSTER
                        genome cluster id
  -m MATRIX, --matrix MATRIX
                        reference SNP matrix
  -t, --tree            build phylogenetic tree

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

## Output
This module generates the following outputs:
* {outdir}/{species_id}.hq_snps: lists of SNPs that passed QC
* {outdir}/{species_id}.ref_freq: reference allele frequency matrix (snps x samples)
* {outdir}/{species_id}.depth: read depth matrix (snps x samples)
* {outdir}/{species_id}.fasta: fasta file containing consensus sequences for species in each sample
* {outdir}/{species_id}.tree: maximum-likelihood phylogenetic tree built from consensus sequences

## Example

## Next steps
[Functionally annotate SNVs] (https://github.com/snayfach/PhyloCNV/blob/master/docs/annotate_snvs.md)

