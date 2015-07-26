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
Download a PhyloCNV database: 
* http://lighthouse.ucsf.edu/phylocnv/phylo_db.tar.gz (140G compressed, 230G expanded)  
  * All 5,952 genome-clusters (31,007 total genomes)  
* http://lighthouse.ucsf.edu/phylocnv/phylo_db.human_gut.tar.gz (7.7G compressed, 11G expanded)  
  * 188 prevalent genome-clusters from the human gut (3,742 total genomes)  
* For more info, see: http://lighthouse.ucsf.edu/phylocnv  
  
And unpack the archive: `tar -zxvf phylo_db.tar.gz`  

### Installation

Download the latest version of the software: https://github.com/snayfach/PhyloCNV/archive/v0.0.1.tar.gz 

Unpack the project: `tar -zxvf PhyloCNV-0.0.1.tar.gz`

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

#### run_phylo_species.py
Use this if you just want to perform metagenomic species (i.e. genome-cluster) profiling
```
usage: run_phylo_species.py [options]

optional arguments:
  -h, --help      show this help message and exit
  -v

Input/Output (required):
  -i INPATH       path to input metagenome in FASTQ/FASTA format. gzip (.gz)
                  and bzip (.bz2) compression supported
  -o OUTBASE      basename for output files: {basename}.abundance,
                  {basename}.summary
  -t TEMP_DIR     path to directory to store temp files (/tmp)

Pipeline Speed (optional):
  -n NREADS       number of reads to use from input metagenome (use all)
  -p THREADS      number of threads to use for database search (1)
  -m              use MicrobeCensus to normalize counts. increases runtime by
                  <=30 additional minutes (False)

Quality control (optional):
  -q MIN_QUALITY  keep reads with quality >= MIN_QUALITY (0)
  -l MIN_LENGTH   keep reads with length >= MIN_LENGTH (0)
  -u MAX_N        keep reads with fraction unknown bases <= MAX_N (1.0)
```

#### run_phylo_cnv.py 
This is the integrated pipeline for estimating the abundance, gene content, and phylogeny of microbial species from a metagenome
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
                        alignment speed/sensitivity (sensitive)
  --reads_align READS_ALIGN
                        # reads to use for pangenome alignment (All)
  --reads_batch RD_BATCH
                        Batch size in # reads. Smaller batch sizes requires
                        less memory, but can take longer to run (5000000)
  --map_pid PID         Minimum percent identity between read and reference
                        (93.0)

SNP detection:
  --snps_mapq SNPS_MAPQ
                        Minimum map quality (20)
  --snps_baseq SNPS_BASEQ
                        Minimum base quality (20)  
```

### Output
Final outputs:  
* genome_clusters.abundance: abundance of 5,952 genome-clusters
* genome_clusters.summary: abundance summary
* coverage: coverage of pan genes within each genome-cluster  
* snps: coverage, allele frequencies, and consensus alleles   

Intermediate files:  
* bam: alignments from mapping metagenome (in batches of reads) to each genome-cluster  
* reassigned: best alignments (>= map_pid) for each read across genome-clusters  
* fastq: reads extracted from reassigned bam files  
* bam_rep: alignments on to representative genome for each genome-cluster  

### Example
Run PhyloCNV using the test FASTQ file:
`run_phylo_cnv.py -1 phylo_cnv/example/SRR413772_1.fastq.gz -D genome_clusters -o phylo_cnv/example --all --verbose`

Run PhyloSpecies using the test FASTQ file: 
`run_phylo_species.py -i phylo_species/example/example_1.fastq.gz -o phylo_species/example/example_1.out --verbose`
