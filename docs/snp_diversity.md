## Population diversity inference

This script will allow you to quantify the genomic diversity of a bacterial population. Diversity can be computed genome-wide, for different classes of genomic sites, or for individual genes. Additionally, diversity can be computed for populations from individual metagenomic samples, or based on data pooled across samples. Finally, there are many options to filter genomic sites.

Before running these scripts, you'll need to have run `merge_midas.py snps` [read more] (https://github.com/snayfach/MIDAS/blob/master/docs/merge_snvs.md).


####USAGE:

`snp_diversity.py --indir /path/to/snps/species_id --out /path/to/output [options]`

####INPUT/OUTPUT OPTIONS:

<b>--indir PATH</b>  
Path to output from 'merge_midas.py snps' for one species.  
Directory should be named according to a species_id and contains files 'snps_*.txt'

<b>--out PATH</b>  
Path to output file  

####DIVERSITY OPTIONS:  
<b>--genomic_type {'genome-wide', 'per-gene'}</b>  
Compute diversity for individual genes or genome-wide (genome-wide)  
	
<b>--sample_type {'pooled-samples', 'per-sample'}</b>  
Compute diversity for individual samples or for pooled data across samples (per-sample)  

<b>--weight_by_depth</b>  
If --sample_type=pooled-samples, weight samples by their mean read depth for the species  

<b>--site_type {'ALL','NC','CDS','1D','2D','3D','4D'}</b>  
Compute diversity using subset of genomic sites sites (ALL)  
ALL=all-sites, NC=non-coding, CDS=coding, XD=X-fold-degenerate-sites  


####SAMPLE SELECTION OPTIONS:  
<b>--sample_depth FLOAT</b>  
Minimum average read depth per sample (0.0)  

<b>--fract_cov FLOAT</b>  
Minimum proportion of genomic sites covered by at least 1 read (0.0)  

<b>--rand_samples INT</b>  
Randomly select INT samples that pass quality control (None)

<b>--keep_samples STR</b>  
Comma-separated list of sample identifiers to use for computing diversity metrics (None)  
Samples are still be subject to other filters  

<b>--exclude_samples STR</b>  
Comma-separated list of sample identifiers to exclude from computing diversity metrics (None)   
Samples are still be subject to other filters  

<b>--max_sample INT</b>  
Maximum number of samples to process (use all)  
Useful for quick tests  


####GENOMIC SITE SELECTION OPTIONS:  
<b>--site_depth INT</b>  
Minimum number of mapped reads per site per sample (2)  

<b>--site_prev FLOAT</b>  
Minimum proportion of samples with least --site_depth at site (0.0)  

<b>--site_maf FLOAT</b>  
Minimum minor allele frequency of site across samples (0.0)  

<b>--site_ratio FLOAT</b>  
Maximum ratio of site-depth to the mean-genome-depth (Inf)  

<b>--site_freq FLOAT</b>  
Minimum combined frequency of reference allele and major alternate allele across samples (0.0)  
Most mapped reads should match the reference allele or major alternate allele.  
Set this value to exclude sites with multiple alternate alleles.  
For example, --site_freq=0.95 excludes these sites:  
  A(ref)=0.80, T(alt1)=0.10, C(alt2)=0.05, G(alt3)=0.05 (Freq_A + Freq_T = 0.90 < 0.95)  
  A(ref)=0.00, C(alt1)=0.80, G(alt2)=0.20, T(alt3)=0.00 (Freq_A + Freq_C = 0.80 < 0.95)  

<b>--max_sites INT</b>  
Maximum number of genomic sites to process (use all)  
Useful for quick tests  

<b>--rand_sites FLOAT</b>  
Randomly select FLOAT proportion of genomic sites that pass quality control  

####READ SELECTION OPTIONS:  
<b>--rand_reads INT</b>  
Randomly select INT reads from each sample for each genomic site  

<b>--replace_reads</b>  
If using --rand_reads, reads are randomly sampled with replacement  


