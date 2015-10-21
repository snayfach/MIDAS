# Tutorial


### Download and installation
PhyloCNV is written in Python and runs on Linux and OSX.  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/install.md)

### Run PhyloCNV
1. Profile species relative abundances in each metagenome  
`run_phylo_cnv species -1 example/fastq_1.gz -o fastq_1.abundances`  
`run_phylo_cnv species -1 example/fastq_2.gz -o fastq_1.abundances`  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/species.md)

2. Profile strain-level gene-content of species in each metagenome  
`run_phylo_cnv genes -1 example/fastq_1.gz -o fastq_1.abundances`  
`run_phylo_cnv genes -1 example/fastq_2.gz -o fastq_1.abundances`  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/cnvs.md)

3. Profile strain-level nucleotide variants of species in each metagenome  
`run_phylo_cnv snvs -1 example/fastq_1.gz -o fastq_1.abundances`  
`run_phylo_cnv snvs -1 example/fastq_2.gz -o fastq_1.abundances`  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/snvs.md)

4. Merge the CNV results for each species across samples
`merge_genes ` 
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_cnvs.md)

5. Merge the SNV results for each species across samples
`merge_genes ` 
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_snvs.md)
