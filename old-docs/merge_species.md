## CAUTION:  This verbatim copy of the original [MIDAS tool](https://github.com/snayfach/MIDAS) documentation is OUT OF DATE for the MIDAS-IGGdb update, and only included as a temporary reference during construction of MIDAS-IGGdb.



## Overview
Merge species abundance files across samples
As input, the script takes a list of sample directories. As output, matrix files are produced with relative abundance matrix, marker gene read-depth, counts of reads mapped to marker genes, and and table of species prevalence.

## Usage
```
Usage: merge_midas.py species <outdir> [options]

positional arguments:
  outdir                Directory for output files

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT              Input to sample directories output by run_midas.py; see '-t' for details
  -t INPUT_TYPE         Specify one of the following:
                          list: -i is a comma-separated list (ex: /samples/sample_1,/samples/sample_2)
                          dir: -i is a directory containing all samples (ex: /samples)
                          file: -i is a file of paths to samples (ex: /sample_paths.txt)
  -d DB                 Path to reference database
                        By default the MIDAS_DB environmental variable is used
  --sample_depth FLOAT  Minimum per-sample marker-gene-depth for estimating species prevalence (1.0)
  --max_samples INT     Maximum number of samples to process.
                        Useful for testing (use all)
```

## Examples

1) provide list of paths to sample directories:  
`merge_midas.py species /path/to/outdir -i /path/to/samples/sample_1,/path/to/samples/sample_2 -t list`  

2) provide directory containing all samples:  
`merge_midas.py species /path/to/outdir -i /path/to/samples -t dir`  

3) provide file containing paths to sample directories:  
`merge_midas.py species /path/to/outdir -i /path/to/samples/sample_paths.txt -t file`  

4) run a quick test:
`merge_midas.py species /path/to/outdir -i /path/to/samples -t dir --max_samples 2`

## Output files

<b>relative\_abundance.txt</b>: relative abundance matrix (columns are samples, rows are species).  
<b>count\_reads.txt</b>: read count matrix (columns are samples, rows are species).  
<b>coverage.txt</b>: genome coverage matrix (columns are samples, rows are species).  
<b>species_prevalence.txt</b>: summary statistics for each species across samples

## Output formats

<b>species_prevalence</b>

* species_id: species identifier      
* species_name: unique species name    
* mean_coverage: average read-depth across samples   
* median_coverage: median read-depth across samples 
* mean_abundance: average relative abundance across samples  
* median_abundance: median relative abundance across samples        
* prevalence: number of samples with >= `MIN_COV`

## Memory usage
* This step takes an insignificant amount of memory  

## Next steps
[Profile strain-level gene content of species](cnvs.md)  
[Profile strain-level SNPs of species](snvs.md)

