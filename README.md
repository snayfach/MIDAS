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
<img src="https://github.com/snayfach/PhyloCNV/blob/master/images/pipeline.pdf" width="700" height="500"/>

