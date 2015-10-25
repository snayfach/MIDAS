## Annotate CNVs
Functionally annotate gene copy number matrix.

* Transform a pangenome matrix (rows are genes) to a function matrix (rows are functions)
* Functions can be one of: [KEGG Pathways] (http://www.genome.jp/kegg/pathway.html), [FIGFam gene families] (http://www.nmpdr.org/FIG/wiki/view.cgi/FIG/FigFam), [Gene Ontology terms] (http://geneontology.org/), or reactions from the [Enzyme Commision database] (http://enzyme.expasy.org/)
* The abundances of genes (presence/absence or copy number) are summed by function_id 
* Genes that have no annotated function are dropped

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
Functionally annotate presence/absence matrix with KEGG Pathways:  
`annotate_snps.py -i 57955.presabs -o 57955.kegg -f kegg`

Functionally annotate copy-number matrix with FIGfams:  
`annotate_snps.py -i 57955.copy_num -o 57955.figfams -f figfams`

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
