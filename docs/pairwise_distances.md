## Build pairwise distance matrix
Compute inter-sample strain-level distances based on a phylogenetic tree, gene presence/absence matrix, or allele frequency matrix.

## Usage
```
usage: pairwise_distances.py [options]

optional arguments:
  -h, --help  show this help message and exit
  -i IN       Input matrix or tree. Accepted file types: .presabs, .copynum,
              .ref_freq, .tree)
  -o OUT      Output distance matrix
```

## Example
Run using defaults:  
`annotate_snps.py -i 57955.presabs -o 57955.kegg -f kegg`

## Output
A matrix where row names and column names are either sample_ids or genome_ids.
Matrix values are distances.


|             | sample_1 | sample_2 | ...  | sample_n |
| :----------:|:-------: | :-------:| :--: | :-------:|
| sample_1    | 0.0      | 0.1      | ...  | 0.4      |
| sample_2    | 0.1      | 0.0      | ...  | 0.2      |
| ...         | ...      | ...      | ...  | ...      |
| sample_n    | 0.4      | 0.2      | ...  | 0.0      |
