## Overview
Merge species abundance across samples

## Usage
```
usage: merge_species.py [options]

Merge species abundance files across samples. Outputs include: a relative
abundance matrix, a genome-coverage matrix, and a table summarizing species
prevalence and abundance across samples

optional arguments:
  -h, --help          show this help message and exit
  -v, --verbose       verbose
  -i INPUT            input to results from 'run_phylo_cnv.py species'. see
                      <intype> for details
  -t {dir,file,list}  input type. 'dir': directory containing species
                      abundance files. example:
                      <directory>/<sample_id>.species 'file': file containing
                      paths to species abundance files. each line in the file
                      should contain the full path to the results for a
                      sample_id. 'list': comma-separated list of paths to
                      phylo_cnv results.
  -o OUTBASE          basename for output files: <outbase>.species_abundance,
                      <outbase>.species_coverage, <outbase>.species_prevalence
  -m MIN_COV          minimum genome-coverage for estimating species
                      prevalence (1.0)                     
```

## Examples

Run using default parameters:  
`merge_species.py -i species -t dir -o outbase`

* `species` is a directory that contains results from `run_phylo_cnv species`
* the name of each file should correspond to a sample_id: \<sample_id>.species

Provide list of file paths:
`merge_species.py -i species/sample_1.species,species/sample_2.species -t list -o outbase`

## Outputs
This script generates three output files:

* **\<out_base>.species_abundance**: relative abundance matrix (columns are samples, rows are species)
* **\<out_base>.species_coverage**: genome coverage matrix (columns are samples, rows are species)
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
[Profile strain-level CNVs of species] (https://github.com/snayfach/PhyloCNV/blob/master/docs/cnvs.md)
[Profile strain-level SNVs of species] (https://github.com/snayfach/PhyloCNV/blob/master/docs/snvs.md)

