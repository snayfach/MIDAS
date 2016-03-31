## Overview
Description: merge results from pan-genome profiling across samples  
Input: list of sample directories  
Output: pan-genome copy-number matrix, presence/absence matrix, and read-depth matrix
        matrixes also created for KEGG, FIGfams, Gene Ontology, and Enzyme Comission (E.C.)

## Usage
```
Usage: merge_genes.py [options]

optional arguments:
  -h, --help            show this help message and exit

Input/Output:
  -i INPUT              input to sample directories output by run_phylo_cnv.py genes
                        see '-t' for details
  -t {list,file,dir}    'list': -i is a comma-separated list of paths to sample directories (ex: /sample1,/sample2)
                        'dir': -i is a  directory containing all samples (ex: /samples_dir)
                        'file': -i is a file containing paths to sample directories (ex: sample_paths.txt)
  -s SPECIES_ID         species identifier
                        a list of prevalent species can be obtained by running 'merge_species.py'
                        a map of species ids to species names can be found in 'ref_db/annotations.txt'
  -o OUTDIR             output directory
  --no_functions        do not write function matrices to OUTDIR

Sample filters (select subset of samples from INPUT):
  --marker_coverage FLOAT
                        minimum coverage per sample across 15 phylogenetic marker genes (1.0)
  --gene_coverage FLOAT
                        minimum coverage per sample across all genes with non-zero coverage (1.0)
  --max_samples INT     maximum number of samples to process. useful for testing (use all)

Presence/Absence:
  --min_copy FLOAT      genes >= MIN_COPY are classified as present
                        genes < MIN_COPY are classified as absent (0.35)
  --cluster_pid {75,80,85,90,95,99}
                        gene family percent identity
                        small values: fewer, larger gene families
                        large values: more, smaller gene families (95)
 ```

## Examples

Examples:
1) Merge results for species 57955. Provide list of paths to sample directories:  
`merge_genes.py -s 57955 -o outdir/57955 -i sample_1,sample_2 -t list`

2) Build matrix for pan-genome genes at lower percent id threshold:  
`merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --cluster_pid 85`

3) Exclude low-coverage samples in output matrix:  
`merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --marker_coverage 5.0`

4) Use lenient threshold for determining gene presence-absence:  
`merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --min_copy 0.1`

5) Just write pan-genome matrices; do not write results for KEGG, FIGfams, GO, or EC:  
`merge_genes.py -s 57955 -o outdir/57955 -i /path/to/samples -t dir --no_functions`


## Outputs
This module generates the following output files:

* **{species_id}.depth**: coverage (i.e. read depth) of each reference gene
* **{species_id}.copynum**: gene copy-number matrix 
  * gene-coverages normalized by the coverage of universal-single-copy genes
* **{species_id}.presabs**: gene presence absence matrix 
  * genes with copy-number >= `min_copy` are called as present and genes below `min_copy` are called as absent (0)


Example gene presence/absence matrix:

| gene_id | sample_1 | sample_2 | ...  | sample_n | genome_1 | ...  | genome_n |
| :----------:|:-------: | :-------:| :--: | :-------:| :-------:| :--: | :-------:|
| 1235786.3.peg.1010       | 1      | 1      | ...  | 1      | 1      | ...  | 1      |
| ...         | ...      | ...      | ...  | ...      | ...      | ...  | ...      |
| 997891.3.rna.48       | 0      | 0      | ...  | 0      | 1      | ...  | 0      |

* Matrix will include reference genomes if `--add_ref` is used
* **gene_id**: gene identifer; `peg` and `rna` indicate coding & non-coding genes respectively

## Memory usage  
* This step takes an insignificant amount of memory


