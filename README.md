# PhyloCNV

####Phylogenetic and gene-copy-number variation of microbial populations

PhyloCNV is an integrated pipeline that estimates bacterial species abundance and strain-level genomic variation from shotgun metagnomes.  Also, population genetic and phylogenetic inferences are performed. PhyloCNV leverage a comprehensive database of ~6,000 bacterial species & ~30,000 genomes.   

## Applications
1. **Profile bacterial species abundance**: rapidly estimate the abundance of ~6,000 bacterial species
2. **Strain-level pan-genome profiling**: identify gene-copy-number variants (CNVs) of strains based on mapping to reference genomes
3. **Single-nucleotide-variant prediction**: identify single-nucleotide variants (SNVs) of strains based on mapping to reference genomes
4. **Phylogenetic inference**: reconstruct the phylogeny of strains from metagenomes and reference genomes
5. **Population genetic inference**: quantify strain-level diversity, differentiation, and selection within and between metagenomes


## Table of Contents
1. [Requirements and Dependencies] (https://github.com/snayfach/PhyloCNV/blob/master/docs/requires.md)
2. [Download, Installation, and Testing] (https://github.com/snayfach/PhyloCNV/blob/master/docs/install.md)
3. [Reference database] (https://github.com/snayfach/PhyloCNV/blob/master/docs/ref_db.md)
4. [Tutorial] (https://github.com/snayfach/PhyloCNV/blob/master/docs/tutorial.md)
5. Scripts to run PhyloCNV on a single sample:
 * [Estimate species abundance] (https://github.com/snayfach/PhyloCNV/blob/master/docs/species.md)
 * [Call gene copy-number variants] (https://github.com/snayfach/PhyloCNV/blob/master/docs/cnvs.md)
 * [Call single-nucleotide variants] (https://github.com/snayfach/PhyloCNV/blob/master/docs/snvs.md)
6. Scripts to merge results across samples:
 * [Merge copy-number variants] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_cnvs.md)
 * [Merge single-nucleotide variants] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_snvs.md)
7. Functionally annotate results:
 * [Annotate copy-number variants] (https://github.com/snayfach/PhyloCNV/blob/master/docs/annotate_cnvs.md)
 * [Annotate single-nucleotide variants] (https://github.com/snayfach/PhyloCNV/blob/master/docs/annotate_snvs.md)
8. [Compare phylogeny, gene-content, and SNPs across samples] (https://github.com/snayfach/PhyloCNV/blob/master/docs/pairwise_distances.md)
9. [Calculate population genetics statistics] (https://github.com/snayfach/PhyloCNV/blob/master/docs/popgen.md)

## Citation
If you use this tool, please cite:
Nayfach, S. and Pollard, KS. PhyloCNV: an integrated, high-resolution pipeline for quantifying strain-level genomic variation from shotgun metagenomes (In Preparation).

## Pipeline
<img src="https://github.com/snayfach/PhyloCNV/blob/master/images/pipeline.jpg" width="600" align="middle"/>   
**An integrated pipeline to estimate bacterial species abundance and strain-level genomic variation from shotgun metagnomes** 
<sub>**A) Metagenome species profiling.** Reads from a metagenomic sample are aligned against a database of phylogenetic marker genes and are assigned to species groups. Mapped reads are used to estimate the genome-coverage and relative abundance of 5,952 genome-clusters. **B) Metagenome pan-genome profiling.** A pan-genome database is dynamically constructed based on the subset of species that are present at high coverage (e.g. >1x) in the metagenome. Reads are mapped to the gene database using Bowtie2. Mapped reads are used to infer gene copy number and gene presence/absence. **C) Single-nucleotide variant prediction.** A representative genome database is constructed, as described in (B). Reads are globally aligned to the genome database using Bowtie2. Mapped reads are used to identify variants, predict consensus alleles, and estimate allele frequencies. **D) Merge results.** For each species, results are merged across one or more samples to generate several outputs, including: a gene presence/absence matrix, an allele frequency matrix, an approximate maximum-likelihood phylogenetic tree.</sub>

## Examples
<img src="https://github.com/snayfach/PhyloCNV/blob/master/images/enrichment.jpg" width="600" align="middle"/>  
**Comparative genomics of *Bacteroides ovatus* strains across host microbiomes**  
<sub> **A)** Presence or absence of genes in the *Bacteroides ovatus* pangenome across human faecal metagenomes. Column colors indicate whether a gene is core (blue; occurs in >95% of samples), auxiliary (red; occurs in 1-95% of samples ), or absent (green; occurs in < 1% of samples). **B)**Gene set enrichment analysis identifies functions overrepresented in the core genome, auxiliary genome, and genes that only occur in reference genomes.</sub>
