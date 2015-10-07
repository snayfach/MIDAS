## Overview
Merge CNVs across samples for a given species

## Usage
```
usage: merge_genes.py [options]

optional arguments:
  -h, --help            show this help message and exit

Input/Output:
  -i IN, --indir IN     input directory
  -o OUT, --outdir OUT  output directory
  -g GENOME_CLUSTER, --genome_cluster GENOME_CLUSTER
                        genome cluster id
  -D DB, --db DB        directory to genome-clusters database

Sample filters:
  --sample_list SAMPLE_LIST
                        file of sample ids to include; each line should
                        contain one id
  --marker_coverage MARKER_COVERAGE
                        min read depth per sample across 15 phylogenetic
                        marker genes (1.0)
  --gene_coverage GENE_COVERAGE
                        min read depth per sample across all genes with non-
                        zero coverage (0.0)
  --max_samples MAX_SAMPLES
                        maximum number of samples to process; useful for
                        testing (use all)

Presence/Absence:
  --min_copy MIN_COPY   genes >= MIN_COPY: classified as present in sample
                        genes < MIN_COPY: classified as absent in sample
                        (0.35)
  --cluster_pid {75,80,85,90,95,99}
                        Gene family percent identity. Small values => fewer,
                        larger gene families. Large values => more, smaller
                        gene families
```

## Output
This module generates the following outputs:
* {outdir}/{species_id}.presabs: gene presence absence matrix (genes x samples)
* {outdir}/{species_id}.copy_num: gene copy number matrix (genes x samples)
* {outdir}/{species_id}.depth: gene depth matrix (genes x samples)

## Example

