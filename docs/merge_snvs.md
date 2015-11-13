## Overview
Merge SNPs across samples for a given species

## Usage
```
usage: merge_snps.py [options]

Merge single-nucleotide variants for an individual species across samples.
Outputs include: a list of high-quality sites, an allele frequency matrix,
consensus sequences for each sample, and a phylogenetic tree

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         verbose

Input/Output:
  -i INPUT              input to results from 'run_phylo_cnv.py snvs'. see
                        <intype> for details
  -t {dir,file,list}    input type. 'dir': directory containing phylo_cnv
                        results. each subdirectory should correspond to a
                        different sample_id. for example:
                        <directory>/<sample_id> 'file': file containing paths
                        to phylo_cnv results. each line in the file should
                        contain the full path to the results for a sample_id.
                        'list': comma-separated list of paths to phylo_cnv
                        results.
  -s SPECIES_ID         species identifier. a list of prevalent species can be
                        obtained by running 'scripts/merge_species.py'. A map
                        of species ids to species names can be found in
                        'ref_db/annotations.txt'
  -o OUTDIR             output directory. output files:
                        <outdir>/<species_id>.hq_snps,
                        <outdir>/<species_id>.ref_freq,
                        <outdir>/<species_id>.depth,
                        <outdir>/<species_id>.fasta,
                        <outdir>/<species_id>.tree

Pipeline options (choose one or more; default=all):
  --snps                identify and store list of hq snps
  --freq                build allele frequency & depth matrixes
  --cons                generate fasta file of consensus sequences
  --tree                build phylogenetic tree

Sample filters
(determine which samples are included in output):
  --sample_depth SAMPLE_DEPTH
                        minimum average read depth per sample (5.0)
  --fract_cov FRACT_COV
                        fraction of reference sites covered by at least 1 read
                        (0.4)
  --max_samples MAX_SAMPLES
                        maximum number of samples to process. useful for quick
                        tests (use all)

Site filters
(determine which reference-genome positions are included in output):
  --site_depth SITE_DEPTH
                        minimum number of mapped reads per site. a high value
                        like 20 will result in accurate allele frequencies,
                        but may discard many sites. a low value like 1 will
                        retain many sites but may not result in accurate
                        allele frequencies (3)
  --site_prev MIN_PREV  site has at least <site_depth> coverage in at least
                        <site_prev> proportion of samples. a value of 1.0 will
                        select sites that have sufficent coverage in all
                        samples. a value of 0.0 will select all sites,
                        including those with low coverage in many samples
                        (0.95)
  --no_fixed            exclude sites with the same consensus allele across
                        samples. this can be useful to reduce the size of the
                        output datasets while retaining most of the
                        information (False)
  --max_sites MAX_SNPS  maximum number of sites to include in output. useful
                        for quick tests (use all)
```

## Examples
Use default parameters and run entire pipeline:  
`merge_snps.py -i snvs -t dir -o Bacteroides_vulgatus -s 57955`

* `snvs` is a directory that contains results from `run_phylo_cnv snvs`
* the name of each subdirectory should correspond to a sample_id
* `57955` is the species identifier for Bacteroides vulgatus

Provide list of sample directories:
`merge_snps.py -i snvs/sample_1,snvs/sample_2 -t list -o Bacteroides_vulgatus -s 57955`

Filter out low-coverage samples (<10x):  
`merge_snps.py -i snvs -t dir -o Bacteroides_vulgatus -s 57955 --sample_depth 10.0`

Exclude variants that have the same consensus allele across all samples:  
`merge_snps.py -i snvs -t dir -o Bacteroides_vulgatus -s 57955 --no_fixed`

Keep variants with missing data in 5% of samples:  
`merge_snps.py -i snvs -t dir -o Bacteroides_vulgatus -s 57955 --snp_prev 0.95`

Identify snps and build allele frequency matrix only:  
`merge_snps.py -i snvs -t dir -o Bacteroides_vulgatus -s 57955 ---snps --freq`

Identify consensus sequences and build phylogenetic tree only:  
`merge_snps.py -i snvs -t dir -o Bacteroides_vulgatus -s 57955 ---cons --tree` 

## Outputs
This module generates the following output files: 

* **{species_id}.hq_snps**: lists of SNPs that passed QC
* **{species_id}.ref_freq**: reference allele frequency matrix (snps x samples)
* **{species_id}.depth**: read depth matrix (snps x samples)
* **{species_id}.fasta**: fasta file containing consensus sequences for species in each sample 
* **{species_id}.tree**: maximum-likelihood phylogenetic tree built from consensus sequences

## Memory usage  
* Memory usage will depend on the number of sites and samples in the resulting output files
* Phylogenetic tree building is the most memory-intensive step. For more info, see: http://www.microbesonline.org/fasttree/#GenomeWide

## Next steps
[Functionally annotate SNVs] (https://github.com/snayfach/PhyloCNV/blob/master/docs/annotate_snvs.md)
[Estimate phylogenetic distance between metagenomes] (https://github.com/snayfach/PhyloCNV/blob/master/docs/pairwise_distances.md)  
[Estimate population genetic parameters]
