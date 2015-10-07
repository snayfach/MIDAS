## Overview
Use Bowtie2 to map metagenomic reads to a database of non-redundant genes from sequenced genomes.
To increase throughput, PhyloCNV builds, and searches reads against, a database that only contains genes from species that are present and abundant in your metagenome.

This module has three main pipeline steps:
* build_db: build bowtie2 database of genes from abundant species
* align: map reads to bowtie2 database
* coverage: compute normalized coverage of genes

## Usage
```
usage: run_phylo_cnv.py genes [options]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose
  --debug               Print out shell commands for debugging purposes
  --remove [{bowtie2_db,bam} [{bowtie2_db,bam} ...]]
                        Remove specified temporary files

Input/Output (required):
  -1 M1                 FASTA/FASTQ file containing 1st mate if paired or
                        unpaired reads
  -2 M2                 FASTA/FASTQ file containing 2nd mate if paired
  -p PROFILE            Path to species profile
  -o OUT                Path to output directory

Pipeline:
  --all                 Run entire pipeline
  --build_db            Build bowtie2 database of pangenomes
  --align               Align reads to pangenome database
  --coverage            Compute coverage of genes in pangenome database

Species to include in pangenome database:
  --gc_topn GC_TOPN     Top N most abundant (None)
  --gc_cov GC_COV       Coverage threshold (None)
  --gc_rbun GC_RBUN     Relative abundance threshold (None)
  --gc_id GC_ID         Identifier of specific genome cluster or comma-
                        separated list of ids (None)

Alignment speed:
  -s {very-fast,fast,sensitive,very-sensitive}
                        Alignment speed/sensitivity (very-sensitive)
  -n READS              # reads to use from input file(s) (use all)
  -t THREADS            Number of threads to use

Computing gene coverage:
  --mapid MAPID         Discard alignments with percent id < MAPID. Higher
                        values indicate fewer mismatches allowed (93.0)
  --aln_cov ALN_COV     Discard alignments where read coverage < ALN_COV.
                        Higher values indicate that reads must be globally
                        covered by alignment (0.70)
```

## Output
* pangenome.bam (intermediate output): alignments of metagenomic reads versus genes
* db (intermediate output): intermediate output that contains fasta and bowtie2 database containing genes
* coverage (final output): directory containing gene coverages for each selected species
  >each file is tab-delimited with a header and three fields: ref_id, raw_coverage, and normalized_coverage
  >normalized coverages are computed by dividing raw_coverages by the median coverage across a panel of 15 universal single-copy genes
  >normalized coverages can be thought of as the average copy number of a gene across a population of cells
* Intermediate outputs can be automatically removed by using the --remove flag

## Example

### Estimate gene coverages only for the most abundant species
run_phylo_cnv.py genes \
-1 PhyloCNV/phylo_cnv/example/example_1.fastq.gz \
-o PhyloCNV/phylo_cnv/example/ex1 \
-p PhyloCNV/phylo_cnv/example/species_abundance.txt \
--all \
--gc_topn 1

### Estimate gene coverages for all species with at least 0.5x coverage
run_phylo_cnv.py genes \
-1 PhyloCNV/phylo_cnv/example/example_1.fastq.gz \
-o PhyloCNV/phylo_cnv/example/ex1 \
-p PhyloCNV/phylo_cnv/example/species_abundance.txt \
--all \
--gc_cov 0.5


## Speed
* Speed, will depend on the number of species you search and the number of reference genomes sequenced per species. 
* For a single species with 1 reference genome, expect ~16,000 reads/second
* Use -n and -t to increase throughput
