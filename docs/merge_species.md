## Overview
Merge species abundance across samples

## Usage
```
Usage: merge_species.py [options]

Description: Merge species abundance files across samples
Input: list of sample directories
Output: relative abundance matrix, genome-coverage matrix, read-count matrix, species prevalence

optional arguments:
  -h, --help          show this help message and exit
  -i INPUT            input to sample directories output by run_phylo_cnv.py species
                      can be a list of directories, a directory containing all samples, or a file with paths
                      see '-t' for details
  -t {list,file,dir}  'list': -i incdicates a comma-separated list of paths to sample directories
                              example: /path/to/samples/sample_1,/path/to/samples/sample_2
                      'dir': -i incdicates a  directory containing all samples
                             example: /path/to/samples
                      'file': -i incdicates a file containing paths to sample directories
                      	   example: /path/to/sample_paths.txt
  -o OUTBASE          basename for output files
  -m MIN_COV          minimum genome-coverage for estimating species prevalence (1.0)

```

## Examples

1) provide list of paths to sample directories:  
`merge_species.py species -i /path/to/samples/sample_1,/path/to/samples/sample_2 -t list -o outbase`

2) provide directory containing all samples:  
`merge_species.py species -i /path/to/samples -t dir -o outbase`

3) provide file containing paths to sample directoriess:  
`merge_species.py species -i /path/to/samples/sample_paths.txt -t file -o outbase`

## Outputs
This script generates three output files:  
* **\<out_base>.relative_abundance**: relative abundance matrix   (columns are samples, rows are species)  
* **\<out_base>.count_reads**: read count matrix (columns are samples, rows are species)  
* **\<out_base>.coverage**: genome coverage matrix (columns are samples, rows are species)  
* **\<out_base>.species_prevalence**: summary statistics for each species across samples  


Example species_prevalence output:

| species_id  | species_name         | mean_coverage | median_coverage  | mean_abundance | median_abundance | prevalence |
| :----------:|:-------:             | :-------:     | :--:             | :-------:      | :-------:        | :--:       |
| 58035       | Bacteroides ovatus   | 2.32          | 2.41             | 0.21           | 0.25             | 2.0        |
| ...         | ...                  | ...           | ...              | ...            | ...              | ...        |
| 54507       | Bacteroides fragilis | 1.44          | 1.44             | 0.12           | 0.12             | 1.0        |

## Memory usage
* This step takes an insignificant amount of memory  

## Next steps
[Profile strain-level gene content of species] (https://github.com/snayfach/PhyloCNV/blob/master/docs/cnvs.md)
[Profile strain-level SNVs of species] (https://github.com/snayfach/PhyloCNV/blob/master/docs/snvs.md)

