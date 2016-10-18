## Strain tracking

These scripts will allow you to identify rare SNPs that discriminate individual strains and to track these SNPs between hosts to elucidate transmission patterns. Before running these scripts, you'll need to have run `merge_midas.py snps` [read more] (https://github.com/snayfach/MIDAS/blob/master/docs/merge_snvs.md).



#### Step 1: identify rare SNPs that disriminate individual strains of a particular species

 * Scan across the entire genome of a patricular species
 * At each genomic site, compute the presence-absence of the four nucleotides across metagenomic samples
 * Identify SNPs (particular nucleotide at a genomic site) that rarely occur across samples
 * Because these SNPs are rare, they serve as good markers of individual strains 

USAGE:

`strain_tracking.py id_markers --indir /path/to/snps/species_id --out species_id.markers [options]`

OPTIONS:

<b>--sample_map PATH </b>  
Path to mapping file that specifes groups of related sample_ids. 
Specify this if there are groups of related samples likely to harbor the same strains (e.g. technical replicates, biological replicates, related individuals).  
The file should be tab-delimited with no header and two columns. Example:  

 sample1 subject1  
 sample2 subject1  
 sample3 subject2  
 sample4 subject2  
 sample5 subject3  
                      
<b>--min_freq FLOAT </b>   
 Minimum allele frequency (proportion of reads) per site for SNP calling (0.10) 
  
<b>--min_reads INT </b>     
Minimum number of reads supporting allele per site for SNP calling (3)  

<b>--max_groups INT </b>   
Discard alleles found in > MAX_GROUPS (1)    
Sample groups are specified by the second field in '--sample_map'  
If '--sample_map' is not specified, each sample_id is treated as a separate group  

<b>--max_sites INT </b>    
Maximum number of genomic sites to process (use all)  
Useful for quick tests

#### Step 2: track rare SNPs between samples and determine transmission 
 * Compute the presence of marker SNPs (identified in Step 1) across metagenomic samples
 * Quantify the number and fraction of marker SNPs that are shared between all pairs of metagenomic samples
 * Based on a sharing cutoff (e.g. 5%), determine if a strain is shared or not

USAGE:

`strain_tracking.py track_markers --indir /path/to/snps/species_id --out species_id.marker_sharing --markers species_id.markers [options]`

OPTIONS:
                      
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
