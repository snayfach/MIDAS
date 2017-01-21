# Metagenomic Intra-Species Diversity Analysis System (MIDAS)


MIDAS is an integrated pipeline that leverages >30,000 reference genomes to estimate bacterial species abundance and strain-level genomic variation, including gene content and SNPs, from shotgun metagnomes. It works particularly well for the human microbiome, or other environments which have many sequenced reference genomes.

## Table of Contents
1. [A brief overview] (docs/overview.md)  
2. [Step-by-step tutorial] (docs/tutorial.md)  
3. [Install or update MIDAS] (docs/install.md)  
4. Reference databse 
 * [Download default database] (docs/ref_db.md)  
 * [Build your own custom database] (docs/build_db.md)
5. Run MIDAS on a single sample:
 * [Estimate species abundance] (docs/species.md)
 * [Predict pan-genome gene content] (docs/cnvs.md)
 * [Call single nucleotide polymorphisms] (docs/snvs.md)
6. Merge MIDAS results across samples:
 * [Merge species abundance] (docs/merge_species.md)  
 * [Merge gene content] (docs/merge_cnvs.md)
 * [Merge SNPs] (docs/merge_snvs.md)
7. Example scripts for analyzing gene content and SNPs:
 * [Strain tracking] (docs/strain_tracking.md)  
 * [Population diversity] (docs/snp_diversity.md)  
 * [Core-genome phylogenetic trees] (docs/snp_trees.md)   
 * [Gene content dynamics] (docs/compare_genes.md)
8. [Citing] (docs/citing.md)