# PhyloCNV
PhyloCNV is an integrated pipeline and for estimating the abundance, gene content, and phylogeny of microbial species from metagnomic data.  PhyloCNV leverage a database of 30,000 Bacterial genomes that have been clustered into species groups using a panel of 30 universal-single-copy genes. 

PhyloCNV consists of three main modules: 
* Species Abundance Estimation: reads are aligned against database of phylogenetic marker genes
* Pan-Gene Presence/Absence & Copy Number: reads are aligned against clusters of reference genomes for abundant species
* Single-Nucleotide-Variant Prediction: reads are aligned against representative genomes for abundant species and samtools mpileup is used to call variants

### Requirements
Python dependencies: 
* Numpy (v1.9.1)
* BioPython (v1.6.2)
* Pysam (v0.8.1)
* MicrobeCensus (v1.0.4)

External packages:
* Bowtie2 (v2.2.4)
* samtools (v1.2)
* bedtools2 (v2.23.0)
* blastn (v2.2.25+)

Tested version numbers are indicated in parenthesis. Other versions may also work.
Linux binaries for external packages are included with this software. Binaries for OSX will be added in the near future.
If included binaries fail to execute, you can compile them on your own machine and place them under: `PhyloCNV/phylo_cnv/bin`

### Reference database
Download PhyloCNV database from: http://lighthouse.ucsf.edu/phylocnv/ 
And unpack the archive: `tar -zxvf phylo_db.tar.gz` 
Note that this archive requires 230G of disk space. Smaller reference databases targeted for specific environments (e.g. marine, human-gut) will be added in the near future

### Installation

Download the latest version of the software: https://github.com/snayfach/PhyloCNV/archive/v0.0.1.tar.gz 

Unpack the project: tar -zxvf PhyloCNV-0.0.1.tar.gz 

Run setup.py. This will install any dependencies: 
`python setup.py install` or
`sudo python setup.py install` to install as a superuser

Alternatively, you can manually install the software.
First, check that required python libraries are installed. You should be able to enter the following command in the python interpreter without getting an error:  
`>>> import Bio.SeqIO`  
`>>> import numpy`
`>>> import pysam`
`>>> import microbe_census`

Next, add the following to your PYTHONPATH environmental variable:  
`export PYTHONPATH=$PYTHONPATH:/path/to/PhyloCNV` or  
`echo -e "\nexport PYTHONPATH=\$PYTHONPATH:/path/to/PhyloCNV" >> ~/.bash_profile` to avoid entering the command in the future

Finally, add the scripts directory to your PATH environmental variable:  
`export PATH=$PATH:/path/to/PhyloCNV/scripts` or  
`echo -e "\nexport PATH=\$PATH:/path/to/PhyloCNV/scripts" >> ~/.bash_profile` to avoid entering the command in the future

Now, you should be able to enter the command into your terminal without getting an error:  
`run_phylo_species.py -h`
`run_phylo_cnv.py -h`

### Usage

```
usage: run_phylo_cnv.py [options]

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -v, --verbose
  -t THREADS, --threads THREADS
                        Number of threads to use

Input/Output (required):
  -1 M1                 FASTQ file containing 1st mate if paired or unpaired
                        reads
  -2 M2                 FASTQ file containing 2nd mate if paired
  -D DB_DIR             Directory of bt2 indexes for genome clusters
  -o OUT                Directory for output files

Pipeline:
  --all                 Run entire pipeline
  --profile             Estimate genome-cluster abundance
  --align               Align reads to genome-clusters
  --map                 Assign reads to mapping locations
  --cov                 Compute coverage of pangenomes
  --extract             Extract mapped reads from bam file & write to FASTQ
  --remap               Re-map reads to representative genomes
  --snps                Run samtools mpileup & estimate SNP frequencies

GC Abundance:
  --reads_gc READS_MS   # reads to use for estimating genome-cluster abundance
                        (5000000)

GC inclusion (choose one):
  --gc_topn GC_TOPN     Top N most abundant (5)
  --gc_cov GC_COV       Coverage threshold (None)
  --gc_rbun GC_RBUN     Relative abundance threshold (None)
  --gc_id GC_ID         Identifier of specific genome cluster (None)
  --gc_list GC_LIST     Comma-separated list of genome cluster ids (None)

Read Alignment/Mapping:
  --align_speed {very-fast,fast,sensitive,very-sensitive}
                        alignment speed/sensitivity (very-sensitive)
  --reads_align READS_ALIGN
                        # reads to use for pangenome alignment (All)
  --reads_batch RD_BATCH
                        Batch size in # reads. Smaller batch sizes requires
                        less memory, but can take longer to run (5000000)
  --map_pid PID         Minimum percent identity between read and reference
                        (93.0)

SNP detection:
  --snps_mapq SNPS_MAPQ
                        Minimum map quality (0)
  --snps_baseq SNPS_BASEQ
                        Minimum base quality (0)  ```
  
### Example
Run PhyloCNV using the test FASTQ file:
`run_phylo_cnv.py -1 phylo_cnv/example/SRR413772_1.fastq.gz -D genome_clusters -o phylo_cnv/example --all --verbose`

Run PhyloSpecies using the test FASTQ file: 
`run_phylo_species.py -i phylo_species/example/example_1.fastq.gz -o phylo_species/example/example_1.out --verbose`
