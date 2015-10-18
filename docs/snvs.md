## Overview
Use Bowtie2 to map metagenomic reads to a database of sequenced genomes.
To increase throughput, PhyloCNV builds, and searches reads against, a database that only contains representative genomes from species that are present and abundant in your metagenome.

This module has three main pipeline steps:
* build_db: build bowtie2 database of genomes from abundant species
* align: map reads to bowtie2 database
* pileup: generate vcf file
* call: parse vcf file and call SNVs for each species

## Usage
```
usage: run_phylo_cnv.py snvs [options]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose
  --debug               Print out shell commands for debugging purposes
  --remove [{bowtie2_db,bam,vcf} [{bowtie2_db,bam,vcf} ...]]
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
  --pileup              Run samtools mpileup
  --call                Call SNPs and format output

Species to include in representative genome database:
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

Read/base filters:
  --mapid MAPID         Discard alignments with percent id < MAPID. Higher
                        values indicate fewer mismatches allowed (93.0)
  --mapq MAPQ           Minimum map quality (20)
  --baseq BASEQ         Minimum base quality (20)
```

## Output
* genomes.bam (intermediate output): alignments of metagenomic reads versus representative genome(s)
* db (intermediate output): contains fasta and bowtie2 database containing genomes
* genomes.vcf (intermediate output): vcf file for all genomes
* snps (final output): directory containing predicted SNVs for each selected species, each file is tab-delimited with a header the fields: 
 * ref_id: scaffold id
 * ref_pos: position on scaffold
 * ref_allele: reference allele
 * alt_allele: alternate allele
 * cons_allele: consensus allele
 * count_alleles: number of alleles observed in metagenome
 * count_ref: count reference alleles observed
 * count_alt: count alternate alleles observed
 * depth: count total reads at ref_pos
 * ref_freq: frequency (0.0 to 1.0) of reference allele
* snps_summary_stats.txt: alignment summary statistics for each species
* Intermediate outputs can be automatically removed by using the --remove flag

## Example

Call SNVs only for the most abundant species
```
run_phylo_cnv.py snvs \
-1 PhyloCNV/phylo_cnv/example/example_1.fastq.gz \
-o PhyloCNV/phylo_cnv/example/ex1 \
-p PhyloCNV/phylo_cnv/example/species_abundance.txt \
--all \
--gc_topn 1
```

Call SNVs for all species with at least 0.5x coverage
```
run_phylo_cnv.py snvs \
-1 PhyloCNV/phylo_cnv/example/example_1.fastq.gz \
-o PhyloCNV/phylo_cnv/example/ex1 \
-p PhyloCNV/phylo_cnv/example/species_abundance.txt \
--all \
--gc_cov 0.5
```
## Speed
* Speed will depend on the number of species you search and the number of sequenced reference genomes per species.
* For a single species with 1 reference genome, expect ~16,000 reads/second
* Use -n and -t to increase throughput
