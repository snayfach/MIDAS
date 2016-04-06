# Tutorial

This step-by-step tutorial will walk you through downloading, installing, and running PhyloCNV

## Download and installation
PhyloCNV is written in Python and runs on Linux and OSX  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/requires.md)  

Download and install the software:  
`git clone https://github.com/snayfach/PhyloCNV`   
`python PhyloCNV/setup.py install`  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/install.md)

Update your PATH:  
`export PATH=$PATH:PhyloCNV/scripts` 

Download the reference database:  
`python PhyloCNV/scripts/download_ref_db.py`   
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/ref_db.md)  

## Run PhyloCNV

More information can be obtained for any of the scripts  using the `-h` flag

Running PhyloCNV can be conceptually broken down in two steps:  
1) run PhyloCNV for each sample: `run_phylo_cnv.py [species, genes, snps]`  
2) merge the results across samples for individual species: `merge_species.py, merge_genes.py, merge_snps.py`

First, move to the example directory and create a new directory to store per-sample output for species, genes, and SNPs:  
`cd PhyloCNV/example`  
`mkdir samples`  

####Profile species abundances

For each sample, run:  
`run_phylo_cnv.py species samples/sample_1 -1 sample_1.fq.gz`  
`run_phylo_cnv.py species samples/sample_2 -1 sample_2.fq.gz`  

* Note that each output directory is named according to a unique sample identifier
* This enables automatically profiling strain-level variation of all species in downstream modules.  
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/species.md)

Merge species abundance across samples:  
`merge_species.py -i samples/sample_1,samples/sample_2 -t list -o example`

* `-i` is a list of sample directories
* `-t list` indicates that `-i` is a list of paths
* `-o` is the basename for output files:
  * \<basename>.count_reads, \<basename>.coverage, \<basename>.relative_abundance: samples x species matrices
  * \<basename>.species_prevalence: summary statistics for species across samples
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_species.md)

####Profile strain-level gene content of abundant species 

For each sample, run:   
`run_phylo_cnv.py genes samples/sample_1 -1 sample_1.fq.gz`  
`run_phylo_cnv.py genes samples/sample_2 -1 sample_2.fq.gz`   

* Note that these directories should already contain species results (ex: samples/sample_1/species)
* Final results are stored in samples/sample_1/genes/coverage
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/cnvs.md)
  
Merge results for individual species across samples:  
`mkdir genes`  
`merge_genes.py -s 57955 -i samples/sample_1,samples/sample_2 -t list -o genes/57955`

* `-s` is a species identifier; 57955 corresponds to Bacteroides vulgatus; by looking at the output from `merge_species.py` we can see that this species was present in sample_1 and sample_2
* `-i` is a list of sample directories
* `-t list` indicates that `-i` is a list of paths
* `-o` is a directory for output files
* The main output, stored in the directory genes/57955, is a pan-genome gene presence/absence matrix (pangenome.presabs)
* PhyloCNV also outputs estimated copy-numbers, and the read depth of all genes (pangenome.copynum, pangenome.depth) 
* These same outputs are included for functional ontologies: KEGG, FIGfams, Gene Ontology, and Enzyme Commission (ex: kegg.copynum)
  *  The values for each functional category are obtained by summing values across pangenome genes that are annotated to functional category
  *  For example a value of 3 for a function in kegg.presabs indicates that there were 3 genes in pangenome.presabs with a value of 1
  * If you just want the pan-genome matrices, use `--no_function`
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_cnvs.md)

####Profile strain-level nucleotide polymorphisms of abundant species  

For each sample, run:  
`run_phylo_cnv.py snps samples/sample_1 -1 sample_1.fq.gz`  
`run_phylo_cnv.py snps samples/sample_2 -1 sample_2.fq.gz`  

* Note that these directories should already contain species results (ex: samples/sample_1/species)
* Final results are stored in samples/sample_1/snps/snps
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/snvs.md)

Merge results for individual species across samples:  
`mkdir snps`  
`merge_snps.py -s 57955 -i samples/sample_1,samples/sample_2 -t list -o snps/57955`

* '57955' is the identifier of *Bacteroides vulgatus*, which was in sample_1 and sample_2   
  * See output from `merge_species.py` for a table listing common species in samples 1 and 2
* The main outputs, stored in the directory 'snps/57955', are: 
  * List of core-genome sites (snps.list), and their annotations (snps.info)
  * Frequency (0.0 to 1.0) of reference alleles for core-genome sites across samples (snps.ref_freq)
  * Core-genome consensus sequences for strains (snps.fasta)
  * Core-genome phylogentic tree for strains (snps.tree)
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_snvs.md)

