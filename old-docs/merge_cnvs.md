## CAUTION:  This verbatim copy of the original [MIDAS tool](https://github.com/snayfach/MIDAS) documentation is OUT OF DATE for the MIDAS-IGGdb update, and only included as a temporary reference during construction of MIDAS-IGGdb.


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
  -i INPUT              Input to sample directories output by run_midas.py
                        Can be a list of directories, a directory containing all samples, or a file with paths
                        See '-t' for details
  -t {list,file,dir}    'list': -i is a comma-separated list of paths to sample directories (ex: /sample1,/sample2)
                        'dir': -i is a  directory containing all samples (ex: /samples_dir)
                        'file': -i is a file containing paths to sample directories (ex: sample_paths.txt)
  -d DB                 Path to reference database
                        By default, the MIDAS_DB environmental variable is used

Species filters (select subset of species from INPUT):
  --min_samples INT     All species with >= MIN_SAMPLES (1)
  --species_id CHAR     Comma-separated list of species ids
  --max_species INT     Maximum number of species to merge. useful for testing (use all)

Sample filters (select subset of samples from INPUT):
  --sample_depth FLOAT  Minimum read-depth across all genes with non-zero coverage (1.0)
  --max_samples INT     Maximum number of samples to process. useful for testing (use all)

Presence/Absence:
  --cluster_pid {75,80,85,90,95,99}
                        In the database, pan-genomes are defined at 6 different % identity clustering cutoffs
                        CLUSTER_PID allows you to quantify gene content for any of these sets of gene clusters
                        By default, gene content is reported for genes clustered at 95% identity (95)
  --min_copy FLOAT      Genes >= MIN_COPY are classified as present
                        Genes < MIN_COPY are classified as absent (0.35)
```

## Examples

Examples:
1) Merge results for all species. Provide list of paths to sample directories:  
`merge_midas.py genes /path/to/outdir -i sample_1,sample_2 -t list`

2) Merge results for one species:  
`merge_midas.py genes /path/to/outdir --species_id Bacteroides_vulgatus_57955 -i sample_1,sample_2 -t list`

3) Quantify genes clusters at 90% identity:  
`merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --cluster_pid 90`

4) Exclude low-coverage samples in output matrix:  
`merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --sample_depth 5.0`

5) Use lenient threshold for determining gene presence-absence:  
`merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --min_copy 0.1`

6) Run a quick test:  
`merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --max_species 1 --max_samples 10`


## Output files
The output of this script generates the following files: 

* <b>genes_copy\_number.txt</b>: gene copy-number matrix (columns are samples, rows are gene families)  
* <b>genes_presence\_absence.txt</b>: gene presence-absence (0/1) matrix (columns are samples, rows are gene families). genes with copy-number >= `MIN_COPY` are called as present and genes below `MIN_COPY` are called as absent
* <b>genes_coverage.txt</b>: gene coverage (i.e. read depth) matrix (columns are samples, rows are gene families)
* <b>genes_summary.txt</b>: alignment summary statistics of genes across samples

## Output formats

<b>genes_summary.txt</b>:

* sample_id: sample identifier      
* pangenome_size: total number of gene families (99% identity clustering cutoff) in reference pangenome 
* covered_genes: number of pangenome gene families with non-zero coverage   
* fraction_covered: fraction of pangenome gene families with non-zero coverage           
* mean_coverage: mean read-depth across gene families with non-zero coverage
* marker_coverage: median read-depth across 15 universal single copy genes

## Memory usage  
* This step takes an insignificant amount of memory


