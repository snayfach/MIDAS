# PhyloCNV
PhyloCNV is an integrated pipeline that estimates bacterial species abundance and strain-level genomic variation from shotgun metagnomes.  
Additionally, PhyloCNV performs population genetic and phylogenetic inferences.  Our software leverage a database of ~30,000 bacterial genomes clustered into ~6,000 species.  PhyloCNV consists of three main modules:
* Species Abundance Estimation
* Pan-Genome Profiling
* Single Nucleotide Variant Prediction

## Table of Contents
1. [Requirements and Dependencies] (https://github.com/snayfach/PhyloCNV/blob/master/docs/requires.md)
2. [Download, Installation, and Testing] (https://github.com/snayfach/PhyloCNV/blob/master/docs/install.md)
3. [Quickstart] (https://github.com/snayfach/PhyloCNV/blob/master/docs/quickstart.md)
4. Scripts to run PhyloCNV on a single sample:
 * [Estimate species abundance] (https://github.com/snayfach/PhyloCNV/blob/master/docs/species.md)
 * [Call gene copy-number variants] (https://github.com/snayfach/PhyloCNV/blob/master/docs/cnvs.md)
 * [Call single-nucleotide variants] (https://github.com/snayfach/PhyloCNV/blob/master/docs/snvs.md)
5. Scripts to merge results across samples:
 * [Merge copy-number variants] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_cnvs.md)
 * [Merge single-nucleotide variants] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_snvs.md)

## Citation
If you use this tool, please cite:
Nayfach, S. and Pollard, KS. PhyloCNV: an integrated, high-resolution pipeline for quantifying strain-level variation from shotgun metagenomes (In Preparation).

## Pipeline
<img src="https://github.com/snayfach/PhyloCNV/blob/master/images/pipeline.jpg"/>
**An integrated pipeline for profiling species abundance and strain-level genomic variation from metagenomes**  
**A)** Metagenome species profiling. Reads from a metagenomic sample are aligned against a database of phylogenetic and are assigned to species groups according to gene-specific and species-level cutoffs. Mapped reads are used to estimate the genome-coverage and relative abundance of 5,952 genome-clusters. **B)** Metagenome pan genome profiling. In the second step of the pipeline, a pan genome database is dynamically constructed based on the subset of genome-clusters that are present at high coverage (e.g. >1x) in the metagenome. Metagenomic reads are locally aligned to the gene database using Bowtie2 and mapped reads are used to obtain pan gene coverages, which are normalized by the median coverage across a panel of 15 universal single-copy genes. **C)** Metagenome single-nucleotide variant prediction. In the third step of the pipeline, a representative genome database is constructed, as described in (B). Metagenomic reads are globally aligned to the genome database using Bowtie2 and mapped reads are used to identify variants, predict consensus alleles, and estimate allele frequencies. **D)** For each genome-cluster, results are merged across one or more samples to generate several outputs. For pan genome analysis, the outputs include a gene presence/absence matrix and a gene copy number matrix. For SNV analysis, the outputs include an allele frequency matrix, core-genome consensus sequences, and an approximate maximum-likelihood phylogenetic tree, which optionally includes phylogenetic placements of sequenced reference genomes.


