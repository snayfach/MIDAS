## Overview
Use HS-BLASTN to map metagenomic reads to a panel of 15 universal, single-copy marker genes (88K total genes). 
These alignments are used to estimate the coverage and relative abundance of 5,952 bacterial species.
Required for downstream steps if selecting species based on their coverage or relative abundance. 
If you already have a list of species ids that you're interested, this step can be skipped.

## Usage
```
usage: run_phylo_cnv.py species [options]

optional arguments:
  -h, --help           show this help message and exit
  -v, --verbose
  -1 M1                FASTA/FASTQ file containing 1st mate if paired or
                       unpaired reads
  -2 M2                FASTA/FASTQ file containing 2nd mate if paired
  -o OUT               Path to output file
  -k                   Keep temporary files, including BLAST output
  -m                   Estimate cellular relative abundance. Requires running
                       MicrobeCensus and takes 20-30 minutes longer to
                       complete.
  -s {fast,sensitive}  Alignment speed/sensitivity (fast)
  -n READS             # reads to use from input file(s) (use all)
  -t THREADS           Number of threads to use
```

## Output
A tab delimited file with a header and three fields:
* cluster_id: species (i.e. genome-cluster) identifier
* coverage: estimated genome-coverage of species in metagenome
* relative_abundance: estimated relative abundance of species in metagenome  
 * If using the option -m, then relative abundances are in terms of the proportion of cells. The total relative abundance across known species will be less than 1.0, and the difference indicates the estimated relative abundance of novel species.  
 * If not using -m, then relative abundances are in terms of the proportion of known species.

## Example
```
run_phylo_cnv.py species \
-1 PhyloCNV/phylo_cnv/example/example_1.fastq.gz \
-o PhyloCNV/phylo_cnv/example/species_abundance.txt
```

## Speed
* ~5,000 reads/second for 100-bp reads when using default parameters
* Using -m will add an additional 20 minutes to runtime
* Use -n and -t to increase throughput
* Note than using -n will result in underestimates of species genome-coverage in the full metagenome
* We found that about 1 million reads was sufficient to precisely estimate species relative abundance for a gut community
