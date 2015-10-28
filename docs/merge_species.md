## Overview
Merge species abundance across samples

## Usage
```
usage: merge_species.py [options]

optional arguments:
  -h, --help   show this help message and exit
  -i IN_DIR    Input directory of species abundance files. Filenames should
               have the format: <in_dir>/<sample_id>.species
  -o OUT_BASE  Basename for output files: <out_base>.species_abundance,
               <out_base>.species_coverage, <out_base>.species_prevalence
  -m MIN_COV   Minimum genome-coverage for estimating species prevalence (1.0)
```

## Examples

Run using default parameters:  
`merge_species.py -i species -o test`

* `species` is a directory that contains results from `run_phylo_cnv species`
* the name of each file should correspond to a sample_id: <sample_id>.species

## Outputs
This script generates three output files:

* **<out_base>.species_abundance**: relative abundance matrix (columns are samples, rows are species)
* **<out_base>.species_coverage**: genome coverage matrix (columns are samples, rows are species)
* **<out_base>.species_prevalence**: summary statistics for each species across samples


Example species_prevalence output:

| species_id  | species_name         | mean_coverage | median_coverage  | mean_abundance | median_abundance | prevalence |
| :----------:|:-------:             | :-------:     | :--:             | :-------:      | :-------:        | :--:       |
| 58035       | Bacteroides ovatus   | 2.32          | 2.41             | 0.21           | 0.25             | 2.0        |
| ...         | ...                  | ...           | ...              | ...            | ...              | ...        |
| 54507       | Bacteroides fragilis | 1.44          | 1.44             | 0.12           | 0.12             | 1.0        |

## Next steps
[Profile strain-level CNVs of species] (https://github.com/snayfach/PhyloCNV/blob/master/docs/cnvs.md)
[Profile strain-level SNVs of species] (https://github.com/snayfach/PhyloCNV/blob/master/docs/snvs.md)

