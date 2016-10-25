## Gene Content Dynamics

These scripts will allow you to compare the gene content of a species between all pairs of metagenomic samples. 

Before running these scripts, you'll need to have run:   
`merge_midas.py genes` [read more] (https://github.com/snayfach/MIDAS/blob/master/docs/merge_snvs.md).


####Command usage:

``` 
gene_sharing.py --indir <PATH> --out <PATH> [options]

  --indir PATH          path to output from 'merge_midas.py genes' for one species
                        directory should be named according to a species_id and contains files 'genes_*.txt')
  --out PATH            path to output file
```

####Options:

```
  --max_genes INT       Maximum number of genes to use. Useful for quick tests (use all)
  --max_samples INT     Maximum number of samples to use. Useful for quick tests (use all)
  --distance {jaccard,euclidean,manhattan}
                        Metric to use for computing distances (jaccard)
  --dtype {presabs,copynum}
                        Data type to use for comparing genes (presabs)
  --cutoff FLOAT        Cutoff to use for determining presence absence (0.35)

```

####Examples:  
1) Run with defaults:  
`gene_sharing.py --indir /path/to/species --out distances.txt`

2) Run a quick test:  
`gene_sharing.py --indir /path/to/species --out distances.txt --max_genes 1000 --max_samples 10`

3) Use a different distance metric:  
`gene_sharing.py --indir /path/to/species --out distances.txt --distance manhattan`

4) Use a lenient cutoff for determining gene presence absence:  
`gene_sharing.py --indir /path/to/species --out distances.txt --cutoff 0.10`

5) Use a strict cutoff for determining gene presence absence:  
`gene_sharing.py --indir /path/to/species --out distances.txt --cutoff 0.75`

####Output format:  
  sample1: first sample identifier  
  sample2: second sample identifier  
  count1: number of present genes in sample1
  count2: number of present genes in sample2  
  count_either: number of genes in sample1 or sample2  
  count_both: number of genes in sample1 and sample2  
  distance: dissimilarity between gene sets 

