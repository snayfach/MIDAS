# Reference database
Contains 31,007 bacterial reference genomes clustered into 5,952 species groups. Species groups are based on 96.5% sequence identity across 30 universal marker genes. These groups correspond to the gold-standard definition of bacterial species based on 95% genome-wide average nucleotide identity (ANI). 

For each of the 5,952 species groups, we identified:

* **Pan-genome**: The set of non-redundant genes (99% identity) across all genomes within species. Used for determining the gene-content of strains in metagenomes
* **Representative genome**: Reference genome that is phylogentically close to other genomes within species. Used for identifying single-nucleotide variants from metagenomes
* **Marker-genes**: 15 universal-single-copy genes from representative genome. Used for rapidly estimating species abundance in metagenome

## Download the reference database
Download the MIDAS reference database of marker genes, pangenomes, and representative genomes:
`python MIDAS/scripts/download_ref_db.py`
This may take several minutes to several hours, depending on your internet speed. The entire database requires ~35 GB of free space to install, once installed it takes about ~17G of free space.  

## A comprehensive genomic resource
<img src="https://github.com/snayfach/MIDAS/blob/master/images/ref_db.jpg" width="350" align="right"/>  **Contruction of MIDAS reference database**
<sub>**A)** 31,007 genomes were hierarchically clustered based on the pairwise identity across a panel of 30 phylogenetic marker genes (pMGs). 5,952 species groups were identified by applying a 96.5% identity cutoff. **B)** Comparison of genome-clusters to annotated species names. Out of 31,007 genomes assigned to a genome-cluster, 5,701 (18%) disagreed with the taxonomic label. Most disagreements are due to genomes lacking annotation at the species level (47%). **C)** Genome-clusters were leveraged to construct three genomic databases to be used for species and strain-level profiling of microbial communities. Arrows denote genes with colors indicating gene families. Non-redundant Pan Genomes: the set of unique (>=99% identity) genes from each genome-cluster. Representative genomes: the most phylogenetically representative genome from each genome-cluster. Phylogenetic marker genes: a set of 15 universal single-copy marker genes from each genome-cluster, which are capable of accurately recruiting metagenomic reads.</sub>  


## Next step
[Run MIDAS on an example dataset] (https://github.com/snayfach/MIDAS/blob/master/docs/tutorial.md)
