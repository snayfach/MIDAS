## CAUTION:  This verbatim copy of the original [MIDAS tool](https://github.com/snayfach/MIDAS) documentation is OUT OF DATE for the MIDAS-IGGdb update, and only included as a temporary reference during construction of MIDAS-IGGdb.


## Population diversity inference

This script will allow you to quantify the genomic diversity of a bacterial population. Diversity can be computed genome-wide, for different classes of genomic sites, or for individual genes. Additionally, diversity can be computed for populations from individual metagenomic samples, or based on data pooled across samples. Finally, there are many options to filter genomic sites.

Before running these scripts, you'll need to have run `merge_midas.py snps` [read more](https://github.com/snayfach/MIDAS/blob/master/docs/merge_snvs.md).


### Usage

```
Usage: snp_diversity.py indir [options]

positional arguments:
  PATH                  path to output from `merge_midas.py snps` for one species
                        directory should be named according to a species_id and contains files 'snps_*.txt')

optional arguments:
  -h, --help            show this help message and exit
  --out PATH            path to output file (/dev/stdout)

Diversity options:
  --genomic_type {genome-wide,per-gene}
                        compute diversity for individual genes or genome-wide (genome-wide)
  --sample_type {per-sample,pooled-samples}
                        compute diversity for individual samples or for pooled reads across samples (per-sample)
  --weight_by_depth     weight data from samples by sequencing depth when --sample_type=pooled-samples
  --rand_reads INT      randomly select N reads from each sample for each genomic site
  --replace_reads       reads drawn with replacement
  --rand_samples INT    randomly select N samples from each genomic site
  --rand_sites FLOAT    randomly select X proportion of high-quality genomic sites
  --snp_maf FLOAT       minor allele frequency cutoff for determining if a site is a SNP (0.01)
  --consensus           call consensus alleles prior to calling SNPs

Sample filters (select subset of samples from INDIR):
  --sample_depth FLOAT  minimum average read depth per sample (0.0)
  --sample_cov FLOAT    fraction of reference sites covered by at least 1 read (0.0)
  --max_samples INT     maximum number of samples to process.
                        useful for quick tests (use all)
  --keep_samples STR    comma-separated list of samples to use for computing diversity metrics.
                        samples will still be subject to other filters
  --exclude_samples STR
                        comma-separated list of samples to exclude from computing diversity metrics.
                        samples will still be subject to other filters

Site filters (select subset of genomic sites from INDIR):
  --site_list PATH      path to file containing newline-delimited list of genomic sites to include.
                        other filters will still apply
  --site_depth INT      minimum number of mapped reads per site (2)
  --site_prev FLOAT     site has at least <site_depth> coverage in at least <site_prev> proportion of samples (0.0)
                        a value of 1.0 will select sites that have sufficent coverage in all samples.
                        a value of 0.0 will select all sites, including those with low coverage in many samples
                        NAs recorded for included sites with less than <site_depth> in a sample
  --site_maf FLOAT      minimum average-minor-allele-frequency of site across samples (0.0)
                        setting this above zero (e.g. 0.01, 0.02, 0.05) will only retain variable sites
                        by default invariant sites are also retained.
  --site_ratio FLOAT    maximum ratio of site-depth to mean-genome-depth (None)
                        a value of 10 will filter genomic sites with 10x greater coverage than the genomic background
  --allele_support FLOAT
                        minimum fraction of reads supporting consensus allele (0.50)
  --locus_type {CDS,RNA,IGR}
                        use genomic sites that intersect: 'CDS': coding genes, 'RNA': rRNA and tRNA genes, 'IGS': intergenic regions
  --site_type {1D,2D,3D,4D}
                        if locus_type == 'CDS', use genomic sites with specified degeneracy: 4D indicates synonymous and 1D non-synonymous sites
  --max_sites INT       maximum number of sites to include in output (use all)
                        useful for quick tests
```
  
### Examples:
1) Quantify within-sample heterogenity genome-wide   
`snp_diversity.py /path/to/snps --genomic_type genome-wide --sample_type per-sample --out /path/to/output`

2) Quantify between-sample heterogenity genome-wide  
`snp_diversity.py /path/to/snps --genomic_type genome-wide --sample_type pooled-sample --out /path/to/output`

3) Quantify between-sample heterogenity per-gene  
`snp_diversity.py /path/to/snps --genomic_type per-gene --sample_type pooled-samples --out /path/to/output`

4) Use downsampling to control the read-depth at each genomic site  
`snp_diversity.py /path/to/snps --genomic_type genome-wide --sample_type per-sample --out /path/to/output`

5) Only quantify diversity at non-synonymous sites  
`snp_diversity.py /path/to/snps --genomic_type genome-wide --sample_type pooled-samples --site_type 1D --locus_type CDS --out /path/to/output`

6) Only quantify diversity at synonymous sites (Note: use the ratio of 5/6 for pN/pS)   
`snp_diversity.py /path/to/snps --genomic_type genome-wide --sample_type pooled-samples --site_type 4D --locus_type CDS --out /path/to/output`
 
7) Quantify SNPs using a different definition of a polymorphism  
`snp_diversity.py /path/to/snps --genomic_type genome-wide --sample_type per-sample --snp_maf 0.05 --out /path/to/output`

8) Run a quick test  
`snp_diversity.py /path/to/snps --max_sites 10000  --out /path/to/output`

