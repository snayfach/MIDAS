## Strain tracking

These scripts will allow you to identify rare SNPs that discriminate individual strains and to track these SNPs between hosts to elucidate transmission patterns.   

Before running these scripts, you'll need to have run:  
`merge_midas.py snps` [read more](https://github.com/snayfach/MIDAS/blob/master/docs/merge_snvs.md).


### Step 1: identify rare SNPs that disriminate individual strains of a particular species

 * Scan across the entire genome of a patricular species
 * At each genomic site, compute the presence-absence of the four nucleotides across metagenomic samples
 * Identify SNPs (particular nucleotide at a genomic site) that rarely occur across samples
 * Because these SNPs are rare, they serve as good markers of individual strains 

#### Command usage:

`strain_tracking.py id_markers --indir <PATH> --out <PATH> [options]`

#### Options:

<b>--samples STR </b>  
Comma-separated list of samples to use for training  
By default, all samples are used
    
<b>--min_freq FLOAT </b>   
 Minimum allele frequency (proportion of reads) per site for SNP calling (0.10)  
  
<b>--min_reads INT </b>     
Minimum number of reads supporting allele per site for SNP calling (3)  

<b>--allele_freq INT </b>   
Maximum occurences of allele across samples (1)   
Setting this to 1 (default) will pick alleles found in exactly 1 sample

<b>--max_sites INT </b>    
Maximum number of genomic sites to process (use all)  
Useful for quick tests
      
#### Examples:
1) Use a subset of sample in SNP matrix for training
`strain_tracking.py id_markers --indir merged_snps/species_id --out species.markers --samples sample1,sample2,sample3`

2) Run a quick test  
`strain_tracking.py id_markers --indir merged_snps/species_id --out species.markers --max_sites 10000`

3) Use strict criteria for pick marker alleles:   
`strain_tracking.py id_markers --indir indir --out outfile --min_freq 0.90 --min_reads 5 --allele_freq 1 `

      
### Step 2: track rare SNPs between samples and determine transmission 
 * Compute the presence of marker SNPs (identified in Step 1) across metagenomic samples
 * Quantify the number and fraction of marker SNPs that are shared between all pairs of metagenomic samples
 * Based on a sharing cutoff (e.g. 5%), determine if a strain is shared or not

#### Command usage:

`strain_tracking.py track_markers --indir /path/to/snps/species_id --out species_id.marker_sharing --markers species_id.markers [options]`

#### Options:
                      
<b>--min_freq FLOAT </b>   
 Minimum allele frequency (proportion of reads) per site for SNP calling (0.10) 
  
<b>--min_reads INT </b>     
Minimum number of reads supporting allele per site for SNP calling (3)  

<b>--max_sites INT </b>    
Maximum number of genomic sites to process (use all)  
Useful for quick tests

<b>--max_samples INT </b>   
Maximum number of samples to process (use all)
Useful for quick tests
