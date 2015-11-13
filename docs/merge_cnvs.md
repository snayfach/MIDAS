## Overview
Merge CNVs across samples for a given species

## Usage
```
usage: merge_genes.py [options]

Merge gene copy-number variants for an individual species across samples.
Outputs include: a gene copy-number matrix, a gene presence/absence matrix,
and a gene read-depth matrix

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         verbose

Input/Output:
  -i INPUT              input to results from 'run_phylo_cnv.py genes'. see
                        <intype> for details
  -t {dir,file,list}    input type. 'dir': directory containing phylo_cnv
                        results. each subdirectory should correspond to a
                        different sample_id. for example:
                        <directory>/<sample_id> 'file': file containing paths
                        to phylo_cnv results. each line in the file should
                        contain the full path to the results for a sample_id.
                        'list': comma-separated list of paths to phylo_cnv
                        results.
  -s SPECIES_ID         species identifier. a list of prevalent species can be
                        obtained by running 'scripts/merge_species.py'. a map
                        of species ids to species names can be found in
                        'ref_db/annotations.txt'
  -o OUTDIR             output directory. output files:
                        <outdir>/<species_id>.copynum,
                        <outdir>/<species_id>.presabs,
                        <outdir>/<species_id>.gene_depth

Sample filters
(determine which samples are included in output):
  --marker_coverage MARKER_COVERAGE
                        minimum read depth per sample across 15 phylogenetic
                        marker genes (3.0)
  --gene_coverage GENE_COVERAGE
                        minimum read depth per sample across all genes with
                        non-zero coverage (0.0)
  --max_samples MAX_SAMPLES
                        maximum number of samples to process. useful for
                        testing (use all)
  --sample_list SAMPLE_LIST
                        file of sample ids to include. each line should
                        contain one id. each id should correspond to a
                        subdirectory under <indir>

Presence/Absence:
  --min_copy MIN_COPY   genes >= MIN_COPY are classified as present and genes
                        < MIN_COPY are classified as absent (0.35)
  --cluster_pid {75,80,85,90,95,99}
                        gene family percent identity. small values: fewer,
                        larger gene families. large values: more, smaller gene
                        families (95)
  --add_ref             include gene presence/absence for reference genomes
 ```

## Examples

Run using default parameters:
`merge_genes.py -i cnvs -t dir -o Bacteroides_vulgatus -s 57955`

* `cnvs` is a directory that contains results from `run_phylo_cnv genes`
* the name of each subdirectory should correspond to a sample_id
* `57955` is the species identifier for Bacteroides vulgatus

Provide list of sample directories:  
`merge_genes.py -i cnvs/sample_1,cnvs/sample_2 -t list -o Bacteroides_vulgatus -s 57955`

Add reference genome pan-genome profiles in output:  
`merge_genes.py -i cnvs -t dir -o Bacteroides_vulgatus -s 57955 --add_ref`

* There will be additional columns in output matrix that correspond to reference genomes
* Enables estimating gene-sharing with reference genomes
* Enables comparison of gene-content variation between metagenomic samples to variation between reference genomes

Adjust the size and diversity of gene-families in output:  
`merge_genes.py -i cnvs -t dir -o Bacteroides_vulgatus -s 57955 --cluster_pid 90`

* Higher values (i.e. 95, 99) will result in additional gene-families in output matrix
* This will result in lower gene-sharing between metagenomes
* Lower values values (i.e. 85, 90) will result in fewer gene-families in output matrix
* This will result in higher gene-sharing between metagenomes

Predict gene presence/absence with higher sensitivity & lower specificity:  
`merge_genes.py -i cnvs -t dir -o Bacteroides_vulgatus -s 57955 --min_copy 0.20`

Predict gene presence/absence with higher specificity & lower sensitivity:  
`merge_genes.py -i cnvs -t dir -o Bacteroides_vulgatuss -s 57955 --min_copy 0.50`

## Outputs
This module generates three output files: 

* **{species_id}.gene_depth**: coverage (i.e. read depth) of each reference gene
* **{species_id}.copy_num**: gene copy-number matrix 
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

## Next steps
[Functionally annotate CNVs] (https://github.com/snayfach/PhyloCNV/blob/master/docs/annotate_cnvs.md)  
[Estimate gene-content distance between metagenomes] (https://github.com/snayfach/PhyloCNV/blob/master/docs/pairwise_distances.md)

