# Tutorial

This step-by-step tutorial will walk you through downloading, installing, and running PhyloCNV

### Download and installation
PhyloCNV is written in Python and runs on Linux and OSX  [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/requires.md)  
You will need to first download the software are reference database [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/install.md)

### Run PhyloCNV
**Profile species abundances in each metagenome**  
`run_phylo_cnv.py species -1 example/sample_1.fq.gz -o example/sample_1.abundances`  
`run_phylo_cnv.py species -1 example/sample_2.fq.gz -o example/sample_2.abundances`  
This enables automatically profiling strain-level variation of all ~6,000 species in downstream modules.  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/species.md)

**Profile strain-level gene content of all abundant species**   
For each sample, run:   
`run_phylo_cnv.py genes -1 sample_1.fq.gz -p sample_1.species -o cnvs/sample_1`  
`run_phylo_cnv.py genes -1 sample_2.fq.gz -p sample_2.species -o cnvs/sample_2`  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/cnvs.md)

For each species, run:  
`merge_genes.py -i cnvs -o B_vulgatus -g 57955`   
57955 is the identifier of *Bacteroides vulgatus*, which was in sample_1 and sample_2  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_cnvs.md)

**Profile strain-level nucleotide variants of abundant species**  
For each sample, run:  
`run_phylo_cnv.py snvs -1 sample_1.fq.gz -p sample_1.species -o snvs/sample_1`  
`run_phylo_cnv.py snvs -1 sample_2.fq.gz -p sample_2.species -o snvs/sample_2`  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/snvs.md)

For each species, run:  
`merge_snps.py -i snvs -o B_vulgatus -g 57955`   
57955 is the identifier of *Bacteroides vulgatus*, which was in sample_1 and sample_2   
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_snvs.md)