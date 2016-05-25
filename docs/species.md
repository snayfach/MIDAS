## Metagenomic species profiling
Use HS-BLASTN to map metagenomic reads to a panel of 15 universal, single-copy marker genes. 
These alignments are used to estimate the coverage and relative abundance of 5,952 bacterial species.
Required for downstream steps if selecting species based on their coverage or relative abundance. 
If you already have a list of species ids that you're interested, this step can be skipped.

## Usage
```
Usage: run_midas.py species outdir [options]

Description:
This script will map metagenomic reads to a database of phylogenetic marker genes using HS-BLASTN
Mapped reads are used estimate the read depth and relative abundance of bacterial species
Reads are mapped according to gene-specific, species-level mapping thresholds (94.5-98% DNA identity)
Reads that map equally well to 2 or more species are probabalistically assigned

positional arguments:
  outdir             Path to directory to store results. Name should correspond to sample identifier.

optional arguments:
  -h, --help         show this help message and exit
  -1 M1              FASTA/FASTQ file containing 1st mate if paired or unpaired reads
  -2 M2              FASTA/FASTQ file containing 2nd mate if paired
  -n MAX_READS       Number of reads to use from input file(s) (use all)
  -t THREADS         Number of threads to use for database search (1)
  --db_type          Reference database. Choices:
                     'phyeco': universal-single-copy protein family database (default)
                     'ssuRNA': 16S ribosomal rna database
  --remove_temp      Remove temporary files, including BLAST output
  --word_size INT    Word size for BLAST search (28)
                     Use word sizes > 16 for greatest efficiency.
  --mapid FLOAT      Discard reads with alignment identity < MAPID
                     By default gene-specific species-level cutoffs are used
                     Values between 0-100 accepted
  --aln_cov FLOAT    Discard reads with alignment coverage < ALN_COV (0.75)
                     Values between 0-1 accepted
  --read_length INT  Trim reads to READ_LENGTH and discard reads with length < READ_LENGTH
                     By default, reads are not trimmed or filtered
```

## Example
1) run with defaults using a paired-end metagenome:  
`run_midas.py species outdir -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz`

2) run using a single-end metagenome with 4 CPUs and only 4M reads:  
`run_midas.py species outdir -1 /path/to/reads_1.fq.gz -t 4 -n 4000000`

3) run with exactly 80 base-pair reads:  
`run_midas.py species outdir -1 /path/to/reads_1.fq.gz --read_length 80`

4) quantify species abundance using a 16S database (not recommended!):  
`run_midas.py species outdir -1 /path/to/reads_1.fq.gz --db_type ssuRNA`


## Output
A tab delimited file with a header and three fields:  
* **species_id**: species (i.e. genome-cluster) identifier  
* **species_name**: consensus genus-species name   
* **count_reads**: number of reads mapped to marker genes  
* **coverage**: estimated genome-coverage of species in metagenome  
* **relative_abundance**: estimated relative abundance of species in metagenome  
        
Example of species abundance table for one sample:

| species_id      | species_name | count_reads | coverage      | relative_abundance  |
| :----------: | :-------------:|:-------------:| :------------------: |
| 60140         | Bacteroides vulgatus | 15000 | 100.67        | 0.45              |
| ...           | ... | ...           |   ...               |
| 60415         | Bacteroides fragilis | 750 |5.05          |   0.09              |


## Memory usage
* < 1.5 Gb for most samples

## Speed
* ~5,000 reads/second for 100-bp reads when using default parameters
* Use `-n` and `-t` to increase throughput
* Note than using `-n` will result in underestimates of species genome-coverage in the full metagenome
* We found that about 1 million reads was sufficient to precisely estimate species relative abundance for a gut community

## Next step
[Merge species abundance across samples] (https://github.com/snayfach/MIDAS/blob/master/docs/merge_species.md)
