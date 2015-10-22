## Metagenomic pan-genome profiling
Use Bowtie2 to map metagenomic reads to a pangenome database, which contains non-redundant genes from sequenced genomes, and estimate the coverage and copy number of these genes.

This module has three main pipeline steps:  

* **build_db**: build bowtie2 pangenome database of non-redundant genes from abundant species  
* **align**: map reads to bowtie2 database  
* **coverage**: compute normalized coverage of genes

## Usage
```
usage: run_phylo_cnv.py genes [options]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose
  --debug               Print out shell commands for debugging purposes
  --remove [{bowtie2_db,bam} [{bowtie2_db,bam} ...]]
                        Remove specified temporary files

Input/Output (required):
  -1 M1                 FASTA/FASTQ file containing 1st mate if paired or
                        unpaired reads
  -2 M2                 FASTA/FASTQ file containing 2nd mate if paired
  -p PROFILE            Path to species profile
  -o OUT                Path to output directory

Pipeline:
  --all                 Run entire pipeline
  --build_db            Build bowtie2 database of pangenomes
  --align               Align reads to pangenome database
  --coverage            Compute coverage of genes in pangenome database

Species to include in pangenome database:
  --gc_topn GC_TOPN     Top N most abundant (None)
  --gc_cov GC_COV       Coverage threshold (None)
  --gc_rbun GC_RBUN     Relative abundance threshold (None)
  --gc_id GC_ID         Identifier of specific genome cluster or comma-
                        separated list of ids (None)

Alignment speed:
  -s {very-fast,fast,sensitive,very-sensitive}
                        Alignment speed/sensitivity (very-sensitive)
  -n READS              # reads to use from input file(s) (use all)
  -t THREADS            Number of threads to use

Computing gene coverage:
  --mapid MAPID         Discard alignments with percent id < MAPID. Higher
                        values indicate fewer mismatches allowed (93.0)
  --aln_cov ALN_COV     Discard alignments where read coverage < ALN_COV.
                        Higher values indicate that reads must be globally
                        covered by alignment (0.70)
```

## Example

Run using defaults:  
```
run_phylo_cnv.py genes -1 sample_1.fq.gz -p sample_1.species -o cnvs/sample_1
```

## Output

The output of this script is a directory with the following files:

* **pangenome.bam** (*intermediate output*): alignments of metagenomic reads versus pangenome database
* **db/** (*intermediate output*): directory that contains local fasta and bowtie2 pangenome 
* **coverage** (*final output*): directory containing gene coverages for each selected species. File names indicate species identifiers.
* intermediate outputs can be removed using `--remove`

Example of gene coverage table for one sample (ex: coverage/57955.cov.gz):

| ref_id      | raw_coverage      | normalized_coverage  |
| :----------: |:-------------:| :------------------: |
| 1235786.3.peg.1010         | 6.78        | 0.78              |
| ...           | ...           |   ...               |
| 997891.3.rna.48         | 3.99          |   0.46              |

Field definitions:  

* **ref_id**: gene identifer; `peg` and `rna` indicate coding & non-coding genes respectively
* **raw_coverage**: (number of aligned base-pairs)/(gene length in base-pairs)
* **normalized_coverage**: (raw_coverage of ref_id)/(median raw_coverage of 15 universal-single-copy genes)
  * this is the estimated average copy number of the gene across the population of cells


## Speed
* Speed will depend on the number of species you search and the number of reference genomes sequenced per species. 
* For a single species with 1 reference genome, expect ~16,000 reads/second
* Use `-n` and `-t` to increase throughput

## Next step
[Merge results across samples] (https://github.com/snayfach/PhyloCNV/blob/master/docs/merge_cnvs.md)

