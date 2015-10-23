## Annotate CNVs
Functionally annotate gene copy number matrix. Aggregate presence/absence values or copy number values by function_id.

## Usage
```
usage: annotate_genes.py [options]

optional arguments:
  -h, --help            show this help message and exit
  -i CNV_MATRIX         Gene CNV matrix. Expected file name:
                        {species_id}.presabs or {species_id}.copynum
  -o FUNCTION_MATRIX    Function CNV matrix
  -f {kegg,figfams,go,ec}
                        kegg=KEGG pathways, figfams=FIGfams, go=Gene Ontology,
                        ec=Enzyme Commission
  -v, --verbose
```

## Example
Run using defaults:  
`annotate_snps.py -i 57955.presabs -o 57955.kegg -f kegg`

## Output
A matrix where row names are function_ids and column names are either sample_ids or genome_ids.
        
Example of function mantrix for one species:

| function_id | sample_1 | sample_2 | ...  | sample_n | genome_1 | ...  | genome_n |
| :----------:|:-------: | :-------:| :--: | :-------:| :-------:| :--: | :-------:|
| 00720       | 1.0      | 1.0      | ...  | 1.0      | 1.0      | ...  | 1.0      |
| ...         | ...      | ...      | ...  | ...      | ...      | ...  | ...      |
| 00920       | 1.0      | 2.0      | ...  | 0.0      | 2.0      | ...  | 3.0      |

Descriptions of function ids can be found in the following files:

* PhyloCNV/ref_db/ontologies/kegg.txt
* PhyloCNV/ref_db/ontologies/figfams.txt
* PhyloCNV/ref_db/ontologies/go.txt
* PhyloCNV/ref_db/ontologies/ec.txt
