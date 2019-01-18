## CAUTION:  This verbatim copy of the original [MIDAS tool](https://github.com/snayfach/MIDAS) documentation is OUT OF DATE for the MIDAS-IGGdb update, and only included as a temporary reference during construction of MIDAS-IGGdb.

# Metagenomic Intra-Species Diversity Analysis System (MIDAS)

MIDAS is an integrated pipeline that leverages >30,000 reference genomes to estimate bacterial species abundance and strain-level genomic variation, including gene content and SNPs, from shotgun metagnomes.

## Table of Contents
1. Getting started
 * [How it works](overview.md)  
 * [Install or update software](install.md)  
 * [Step-by-step tutorial](tutorial.md)  
 * [Docker-image](https://github.com/FredHutch/docker-midas)  
2. Reference databse
 * [Download default database](ref_db.md)  
 * [Build your own custom database](build_db.md)
4. Run MIDAS on a single sample:
 * [Estimate species abundance](species.md)
 * [Predict pan-genome gene content](cnvs.md)
 * [Call single nucleotide polymorphisms](snvs.md)
5. Merge MIDAS results across samples:
 * [Merge species abundance](merge_species.md)  
 * [Merge gene content](merge_cnvs.md)
 * [Merge SNPs](merge_snvs.md)
6. Example scripts for analyzing gene content and SNPs:
 * [Core-genome phylogenetic trees](snp_trees.md)
 * [Population diversity](snp_diversity.md)
 * [Strain tracking](strain_tracking.md)      
 * [Gene content dynamics](compare_genes.md)
7. [Citing](citing.md)
