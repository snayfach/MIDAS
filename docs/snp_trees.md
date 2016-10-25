## Core-Genome Phylogenetic Trees

This script will allow you to build strain-level phylogenetic trees using consensus-alleles found in the core-genome.   

The core-genome of a species is identified directly from the data by looking for genomic sites in the representative genome that have high coverage across multiple metagenomes.  

Before running these scripts, you'll need to have run `merge_midas.py snps` [read more] (https://github.com/snayfach/MIDAS/blob/master/docs/merge_snvs.md).

### Step 1. Call consensus alleles

####Command usage:  

``` 
call_consensus.py --indir <PATH> --out <PATH> [options]

  --indir PATH          path to output from 'merge_midas.py snps' for one species
                        directory should be named according to a species_id and contains files 'snps_*.txt')  
  --out PATH            path to output file
```

####Options for selecting samples:

```
  --sample_depth FLOAT  minimum average read depth per sample (0.0)  
  --sample_cov FLOAT    fraction of reference sites covered by at least 1 read (0.0)  
  --max_samples INT     maximum number of samples to process.  
                        useful for quick tests (use all)  
  --keep_samples STR    comma-separated list of samples to include  
                        samples will still be subject to other filters  
  --exclude_samples STR
                        comma-separated list of samples to exclude.  
                        samples will still be subject to other filters  
```

####Options for selecting genomic sites:

```
  --site_depth INT      minimum number of mapped reads per site (2)
  --site_prev FLOAT     site has at least <site_depth> coverage in at least <site_prev> proportion of samples (0.0)
                        a value of 1.0 will select sites that have sufficent coverage in all samples.
                        a value of 0.0 will select all sites, including those with low coverage in many samples 
                        NAs recorded for included sites with less than <site_depth> in a sample 
  --site_maf FLOAT      minimum average-minor-allele-frequency of site across samples (0.0)
                        setting this above zero (e.g. 0.01, 0.02, 0.05) will only retain variable sites
                        by default invariant sites are also retained.
  --site_ratio FLOAT    maximum ratio of site-depth to mean-genome-depth (None)
                        a value of 10 will filter genomic sites with 10x high coverage than the genomic background
  --max_sites INT       maximum number of sites to include in output (use all)
                        useful for quick tests
```

####Example command:  
```
call_consensus.py --indir /path/to/snps --out consensus.fa --site_maf 0.01 --site_depth 5 --site_prev 0.90 --sample_depth 10 --sample_cov 0.40 --site_ratio 
```

This command will build a multi-FASTA of core-genome sequences   
-core-genome sites defined as >=5 reads in >=90% of samples  
-use only variable positions (>=1% minor allele frequency across samples)  
-only include samples with sufficient data (>=10x mean-depth, >=40% of sites with >=1 mapped read)  
-exclude sites with abnormal depth (>5x mean-depth or <1/5 mean-depth)  



### Step 2. Build phylogenetic tree
Now simply use your favorite tool to build the phylogenetic tree.

Example:  
Download FastTree [here] (http://www.microbesonline.org/fasttree)  
And, run: `FastTree -gtr -nt < consensus.fa > consensus.tree `


### Step 3. Visualize tree
From a web-browser: [iTOL] (http://itol.embl.de/)  
Or, from R: [ape] (https://cran.r-project.org/web/packages/ape)

