## Overview
Description: perform pooled-sample core-genome SNP calling. For each species, this script will pool data across multiple samples. Based on this pooled data, it will identify core genomic sites present in most samples. Next it will determine if any of these core-sites are bi-allelic SNPs (i.e. there are two alleles at >5% frequency). Take a close look at the options; there are different ways of doing this step.


The pipeline can be broken down into the following steps:

  * Take MIDAS output files from multiple samples
  * Identify species to process (based on user criterea, e.g. min # of samples)
  * Scan across the representative genome of each species
  * Pool nucleotide variants from all metagenomic samples 
  * Determine if genomic site is in the core-genome (e.g. non-zero depth in >95% of samples)
  * If yes, identify the major and minor alleles among the pooled reads
  * Determine if there is a SNP present (e.g. minor allele frequency >1%)
  * Annotate genomic site (e.g. CDS/RNA, gene identifier, SYN/NS, etc.)
  * Write core-genome SNPs to matrix files
  
## Usage
```
Usage: merge_midas.py snps outdir [options]

positional arguments:
  outdir                Directory for output files. a subdirectory will be created for each species_id

Input/Output:
  -i INPUT              Input to sample directories output by run_midas.py; see '-t' for details
  -t INPUT_TYPE         Specify one of the following:
                          list: -i is a comma-separated list (ex: /samples/sample_1,/samples/sample_2)
                          dir: -i is a directory containing all samples (ex: /samples)
                          file: -i is a file of paths to samples (ex: /sample_paths.txt)
  -d DB                 Path to reference database
                        By default, the MIDAS_DB environmental variable is used

Presets:
  --core_snps           Same as: --snp_type bi --site_depth 1 --site_ratio 2.0 --site_prev 0.95 (default)
  --core_sites          Same as: --snp_type any --site_depth 1 --site_ratio 2.0 --site_prev 0.95
  --all_snps            Same as: --snp_type bi --site_prev 0.0
  --all_sites           Same as: --snp_type any --site_prev 0.0

Species filters (select subset of species from INPUT):
  --min_samples INT     All species with >= MIN_SAMPLES (1)
  --species_id CHAR     Comma-separated list of species ids
  --max_species INT     Maximum number of species to call SNPs for (all with >= 1 sample)

Sample filters (select subset of samples from INPUT):
  --sample_depth FLOAT  Minimum average read depth per sample (5.0)
  --fract_cov FLOAT     Fraction of reference sites covered by at least 1 read (0.4)
  --max_samples INT     Maximum number of samples to process. useful for quick tests (use all)
  --all_samples         Include all samples in output

Site filters (select subset of genomic sites from INPUT):
  --snp_type  [ ...]    Specify one or more of the following:
                          mono: keep sites with 1 allele > ALLELE_FREQ
                          bi: keep sites with 2 alleles > ALLELE_FREQ (default)
                          tri: keep sites with 3 alleles > ALLELE_FREQ
                          quad: keep sites with 4 alleles > ALLELE_FREQ
                          any: keep sites regardless of observed alleles
  --allele_freq FLOAT   Minimum frequency for calling an allele present (0.01)
                        Values > 0.0 and < 0.5 are accepted.
                        Ex: --snp_type=bi --allele_freq=0.01 keeps bi-allelic SNPs with a minimum frequency of 1%
  --site_depth INT      Minimum number of reads mapped to genomic site (1)
                        Used in combination with --site_prev to determine if site is in core-genome
  --site_ratio FLOAT    Maximum ratio of site depth to genome depth (2.0)
                        This filter helps to eliminate genomic sites with abnormally high read depth
  --site_prev FLOAT     Minimum fraction of sample where genomic site is >= SITE_DEPTH and <= SITE_RATIO (0.95)
                        A high value selects for sites in the core-genome (i.e. present in nearly all strains).
                        A low value includes sites in variable regions and/or with abnormally high read depth
  --max_sites INT       Maximum number of sites to include in output (use all). Useful for quick tests
```
## Examples

1) Call core-genome SNPs all species using defaults:  
`merge_midas.py snps /path/to/outdir -i sample_1,sample_2 -t list`  

2) Call core-genome SNPs for one species:
`merge_midas.py snps /path/to/outdir --species_id Bacteroides_vulgatus_57955 -i sample_1,sample_2 -t list`

3) Include all core-genome sites in output files (useful for quantifying SNP density or comparing SNP levels between species):
merge_midas.py snps /path/to/outdir -i /path/to/samples -t dir --core-sites

4) Include ALL sites in output files, regardless of read-depth or variation:
merge_midas.py snps /path/to/outdir -i /path/to/samples -t dir --core-sites

5) Only use samples with >15x average depth and only use sites covered by >=10 reads in at least >=95% of samples:  
`merge_midas.py snps /path/to/outdir -i /path/to/samples -t dir --sample_depth 15 --site_depth 10 --site_prev 0.95`

6) Run a quick test:  
`merge_midas.py snps /path/to/outdir -i /path/to/samples -t dir --max_species 1 --max_samples 10 --max_sites 1000`

## Output files

* <b>snps_freq.txt:</b> frequency of minor allele per genomic site and per sample; a value of 1.0 indicates that all reads matched the minor allele for site-sample; the major (most common) and minor allele (2nd most common) are determined from pooled reads across ALL samples; see: `snps_info.txt` for details on the major, minor, and reference alleles
* <b>snps_depth.txt:</b> number of reads mapped to genomic site per sample; only accounts for reads matching either major or minor allele
* <b>snps_info.txt:</b>; metadata for genomic site; see below for more information
* <b>snps_summary.txt:</b> alignment summary statistics per sample; see below for more information
* <b>snps_log.txt:</b> log file containing parameters used

## Output formats

<b>snps\_freq.txt and snps\_depth.txt</b>: tab-delimited matrix files; field names are sample ids; row names are genome site ids; see: `snps_info.txt` for details on each genomic site
  
<b>snps_info.txt</b>

  * site_id: incrementing integer field
  * ref_id: identifier of scaffold in representative genome
  * ref_pos: position of site on scaffold
  * ref_allele: allele in reference genome
  * major_allele: most common allele in metagenomes
  * minor_allele: second most common allele in metagenomes
  * count_samples: number of metagenomes where site was found
  * count_a: count of A allele in pooled metagenomes
  * count_c: count of C allele in pooled metagenomes
  * count_g: count of G allele in pooled metagenomes
  * count_t: count of T allele in pooled metagenomes
  * locus_type: CDS (site in coding gene), RNA (site in non-coding gene), IGR (site in intergenic region)
  * gene_id: gene identified if locus type is CDS, or RNA
  * snp_type: indicates the number of alleles observed at site (mono,bi,tri,quad); observed allele are determined by `--snp_maf` flag
  * site_type: indicates degeneracy: 1D, 2D, 3D, 4D
  * amino_acids: amino acids encoded by 4 possible alleles

<b>snps_summary.txt</b>

  * sample_id: sample identifier
  * genome_length: number of base pairs in representative genome
  * covered_bases: number of reference sites with at least 1 mapped read
  * fraction_covered: proportion of reference sites with at least 1 mapped read
  * mean_coverage: average read-depth across reference sites with at least 1 mapped read
  * aligned_reads: number of reads that aligned BEFORE quality filtering
  * mapped_reads: number of reads that aligned AFTER quality filtering

## Memory usage  
* This step takes an insignificant amount of memory
