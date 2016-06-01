## Overview
Description: merge results from pan-genome profiling across samples
Input: list of sample directories
Output: pan-genome copy-number matrix, presence/absence matrix, and read-depth matrix
        matrixes also created for KEGG, FIGfams, Gene Ontology, and Enzyme Comission (E.C.)

## Usage
```
Usage: merge_midas.py genes outdir [options]

positional arguments:
  outdir                directory for output files. a subdirectory will be created for each species_id

optional arguments:
  -h, --help            show this help message and exit

Input/Output:
  -i INPUT              input to sample directories output by run_midas.py
                        see '-t' for details
  -t {list,file,dir}    'list': -i is a comma-separated list of paths to sample directories (ex: /sample1,/sample2)
                        'dir': -i is a  directory containing all samples (ex: /samples_dir)
                        'file': -i is a file containing paths to sample directories (ex: sample_paths.txt)

Species filters (select subset of species from INPUT):
  --min_samples INT     all species with >= MIN_SAMPLES (1)
  --species_id CHAR     comma-separated list of species ids
                        a list of prevalent species can be obtained by running 'merge_midas.py species'
                        a map of species ids to species names can be found in 'ref_db/annotations.txt'
  --max_species INT     maximum number of species to merge. useful for testing (use all)

Sample filters (select subset of samples from INPUT):
  --sample_depth FLOAT  minimum coverage per sample across all genes with non-zero coverage (1.0)
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
1) Merge results for all species. Provide list of paths to sample directories:
`merge_midas.py genes /path/to/outdir -i sample_1,sample_2 -t list`

2) Merge results for one species (id=57955):
`merge_midas.py genes /path/to/outdir --species_id 57955 -i sample_1,sample_2 -t list`

3) Build matrix for pan-genome genes at lower percent id threshold:
`merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --cluster_pid 85`

4) Exclude low-coverage samples in output matrix:
`merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --sample_depth 5.0`

5) Use lenient threshold for determining gene presence-absence:
`merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --min_copy 0.1`

6) Run a quick test:
`merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --max_species 1 --max_samples 10`


## Outputs
The output of this script generates the following files: 

* **genes_copy_number.txt**: gene copy-number matrix (columns are samples, rows are gene families)
* **genes_presence_absence.txt**: gene presence-absence (0/1) matrix (columns are samples, rows are gene families). genes with copy-number >= `MIN_COPY` are called as present and genes below `MIN_COPY` are called as absent
* **genes_coverage.txt**: gene coverage (i.e. read depth) matrix (columns are samples, rows are gene families)
* **genes_info.txt**: detailed information of genes 
* **genes_summary.txt**: alignment summary statistics of genes across samples

**genes_info.txt** output format:

* gene_id: identifier of 99% identity gene family
* family_id: mapping to gene family clustered at CLUSTER_PID
* function_id: identifier of function     
* function_db: database (kegg, figfam, go, ec) corresponding to function_id

**genes_summary.txt** output format:

* sample_id: sample identifier      
* pangenome_size: total number of gene families (99% identity clustering cutoff) in reference pangenome 
* covered_genes: number of pangenome gene families with non-zero coverage   
* fraction_covered: fraction of pangenome gene families with non-zero coverage           
* mean_coverage: mean read-depth across gene families with non-zero coverage
* marker_coverage: median read-depth across 15 universal single copy genes

## Memory usage  
* This step takes an insignificant amount of memory


