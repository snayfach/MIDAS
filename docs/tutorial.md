# Tutorial

This step-by-step tutorial will walk you through downloading, installing, and running PhyloCNV

## Download and installation
PhyloCNV is written in Python and runs on Linux and OSX  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/requires.md)  

Download and install the software:  
`git clone https://github.com/snayfach/PhyloCNV`   
`python PhyloCNV/setup.py install`  
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/install.md)

Download the reference database:  
`python PhyloCNV/scripts/download_ref_db.py`   
[read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/ref_db.md)  

## Run PhyloCNV

The following directions assume that PhyloCNV/scripts has been added to your PATH  

More information can be obtained for any of the scripts  using the `-h` flag

Move to the example directory:  
`cd PhyloCNV/example`  

####Profile species abundances

Create a directory to store species abundance files:  
`mkdir species`

For each sample, run:  
`run_phylo_cnv.py species -1 sample_1.fq.gz -o species/sample_1.species --verbose`
`run_phylo_cnv.py species -1 sample_2.fq.gz -o species/sample_2.species --verbose`  

* This enables automatically profiling strain-level variation of all ~6,000 species in downstream modules.  
* Output file names (-o) should have the format \<sample_id>.species
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/species.md)

Merge results across samples:  
`merge_species.py -i species -o merged`

* `-i` is the directory containing \<sample_id>.species files
* `-o` is the basename for output files:
  * \<basename>.species_prevalence: table mapping species ids to species names & the frequency of each species across samples
  * \<basename>.species_abundance: matrix of samples x species
  * \<basename>.species_coverage: matrix of samples x species
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_species.md)

####Profile strain-level gene content of abundant species 

Create a directory to store output files:  
`mkdir cnvs`

For each sample, run:   
`run_phylo_cnv.py genes -1 sample_1.fq.gz -p sample_1.species -o cnvs/sample_1  --verbose`  
`run_phylo_cnv.py genes -1 sample_2.fq.gz -p sample_2.species -o cnvs/sample_2  --verbose`   

* 'cnvs' is the output directory where we will store our results
* the output directory is created if it doesn't already exist  
* 'cnvs/sample_1' and 'cnvs/sample_2' are subdirectories for each sample that contain the output files
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/cnvs.md)

Merge results for each species across samples:  
`merge_genes.py -i cnvs -o Bacteroides_vulgatus -s 57955`

* this script takes a directory as input, where subdirectories correspond to sample_ids 
* '57955' is the identifier of *Bacteroides vulgatus*, which was in sample_1 and sample_2   
  * See merged.species_prevalence for a table listing common species in samples 1 and 2
* The main output, stored in the directory 'Bacteroides_vulgatus', is a gene presence/absence matrix  
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_cnvs.md)

####Profile strain-level nucleotide variants of abundant species  
For each sample, run:  
`run_phylo_cnv.py snvs -1 sample_1.fq.gz -p sample_1.species -o snvs/sample_1`  
`run_phylo_cnv.py snvs -1 sample_2.fq.gz -p sample_2.species -o snvs/sample_2`  

* 'snvs' is the output directory where we will store our results
* the output directory is created if it doesn't already exist  
* 'snvs/sample_1' and 'snvs/sample_2' are subdirectories for each sample that contain the output files
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/snvs.md)

Merge results for each species across samples:    
`merge_snps.py -i snvs -o Bacteroides_vulgatus -s 57955`   

* '57955' is the identifier of *Bacteroides vulgatus*, which was in sample_1 and sample_2   
  * See merged.species_prevalence for a table listing common species in samples 1 and 2
* The main outputs, stored in the directory 'Bacteroides_vulgatus', are: 
  * SNV frequency matrix
  * Core-genome consensus sequences for strains
  * Core-genome phylogentic tree for strains  
* [read more...] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_snvs.md)


## Next steps
[Functionally annotate CNVs] (https://github.com/snayfach/PhyloCNV/blob/master/docs/annotate_cnvs.md)  
[Functionally annotate SNVs] (https://github.com/snayfach/PhyloCNV/blob/master/docs/annotate_snvs.md)  
[Estimate gene-content distance between metagenomes] (https://github.com/snayfach/PhyloCNV/blob/master/docs/pairwise_distances.md)  
[Estimate phylogenetic distance between metagenomes] (https://github.com/snayfach/PhyloCNV/blob/master/docs/pairwise_distances.md)    

