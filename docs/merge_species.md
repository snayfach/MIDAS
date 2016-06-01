## Overview
Merge species abundance files across samples
Input: list of sample directories
Output: relative abundance matrix, genome-coverage matrix, read-count matrix, species prevalence

## Usage
```
Usage: merge_midas.py species outdir [options]

optional arguments:
  -h, --help          show this help message and exit
  -i INPUT            input to sample directories output by run_midas.py species
                      can be a list of directories, a directory containing all samples, or a file with paths
                      see '-t' for details
  -t {list,file,dir}  'list': -i incdicates a comma-separated list of paths to sample directories
                              example: /path/to/samples/sample_1,/path/to/samples/sample_2
                      'dir': -i incdicates a  directory containing all samples
                             example: /path/to/samples
                      'file': -i incdicates a file containing paths to sample directories
                      	   example: /path/to/sample_paths.txt
  -m MIN_COV          minimum genome-coverage for estimating species prevalence (1.0)
```

## Examples

1) provide list of paths to sample directories:  
`merge_midas.py species /path/to/outdir -i /path/to/samples/sample_1,/path/to/samples/sample_2 -t list`  

2) provide directory containing all samples:  
`merge_midas.py species /path/to/outdir -i /path/to/samples -t dir`  

3) provide file containing paths to sample directoriess:  
`merge_midas.py species /path/to/outdir -i /path/to/samples/sample_paths.txt -t file`  

## Outputs
This script generates the following files:  

* **relative_abundance.txt**: relative abundance matrix (columns are samples, rows are species)
* **count_reads.txt**: read count matrix (columns are samples, rows are species)
* **coverage.txt**: genome coverage matrix (columns are samples, rows are species)
* **species_prevalence.txt**: summary statistics for each species across samples

species_prevalence output format:

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
[Profile strain-level gene content of species] (https://github.com/snayfach/MIDAS/blob/master/docs/cnvs.md)  
[Profile strain-level SNVs of species] (https://github.com/snayfach/MIDAS/blob/master/docs/snvs.md)

