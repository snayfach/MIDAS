## Annotate SNVs
Functionally annotate genomic sites and variants

## Usage
```
usage: annotate_snps.py [options]

optional arguments:
  -h, --help            show this help message and exit
  -i IN, --snp_list IN  SNP list
  -o OUT, --annotations OUT
                        SNP annotations
  -m MAX_SNPS, --max_snps MAX_SNPS
                        Maximum # of SNPs to annotate (all)  
  -v, --verbose
```

## Example
Run using defaults:  
`annotate_snps.py -i 57955.hq_snps -o 57955.snp_annotations`


## Output
A tab delimited file with a header and six fields:  
        
Example of snp annotation table for one species:

| ref_id            | ref_pos | count_alt | alt_allele | site_type | snp_type | gene_id |
| :----------:      |:-------:| :-------: | :-------:  | :-------: |:-------:|:-------:|
| accn\|ASTZ01000001 | 1       | 0         | NA         | NA        | NA      | NA      |
| accn\|ASTZ01000001 | 2       | 1         | NA         | NC        | NA      | NA      |
| accn\|ASTZ01000001 | 3       | 1         | T          | ND        | NS      | NA      |    
| accn\|ASTZ01000001 | 4       | 1         | T          | 2D        | SYN      | NA      |    
| accn\|ASTZ01000001 | 5       | 1         | A         | 2D        | NS      | NA      |     
| accn\|ASTZ01000001 | 6       | 1         | C         | 3D        | SYN      | NA      |
| accn\|ASTZ01000001 | 7       | 1         | C         | 3D        | NS      | NA      |
| accn\|ASTZ01000001 | 8       | 1         | C         | 4D        | SYN      | NA      |

Field definitions:  

* **ref_id**: scaffold id of bacterial genome
* **ref_pos**: genomic position on scaffold id
* **count_alt**: number of alternate alleles observed across samples
* **alt_allele**: alternate allele {A,T,C,G}
* **site_type**: NA (no alternate allele), NC (non-coding), ND (non-degenerate), 2D (2-fold degenerate), 3D (3-fold degenerate), 4D (4-fold degenerate)
* **snp_type**: NA (not a snp), SYN (synonymous change), NS (non-synonymous change)


