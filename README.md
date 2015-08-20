# PhyloCNV
PhyloCNV is an integrated pipeline and for estimating the abundance, gene content, and phylogeny of microbial species from metagnomic data.  PhyloCNV leverage a database of 30,000 Bacterial genomes that have been clustered into species groups using a panel of 30 universal-single-copy genes. 

PhyloCNV consists of three main modules: 
* Species Abundance Estimation  
 -rapidly map reads to db of universal genes & probabalistally assign reads to species groups  
 -estimate genome coverage of species-groups   

* Pan Genome Alignment and Coverage  
 -build a bowtie2 database of pangenomes from abundant species    
 -use bowtie2 to map reads to pangenome database  
 -compute normalized coverage of genes  

* Single Nucleotide Variant Prediction  
-build a bowtie2 database of representative genomes from abundant species  
-use bowtie2 to map reads to genome database  
-call SNVs and estimate allele frequencies     

### Requirements
Python dependencies (installed via setup.py): 
* Numpy (v1.9.1)
* BioPython (v1.6.2)
* Pysam (v0.8.1)

External packages (included):
* Bowtie2 (v2.2.4)
* samtools (v1.2)
* blastn (v2.2.25+)

Tested version numbers are indicated in parenthesis. Other versions may also work.
Linux binaries for external packages are included with this software. If included binaries fail to execute, you can compile them on your own machine and place them under: `PhyloCNV/phylo_cnv/bin/Linux` or `PhyloCNV/phylo_cnv/bin/Darwin`

### Reference database
Download a PhyloCNV database: 
* http://lighthouse.ucsf.edu/phylocnv/genome_clusters.tar.gz (17G expanded)  
* For more info, see: http://lighthouse.ucsf.edu/phylocnv  
  
And unpack the archive: `tar -zxvf genome_clusters.tar.gz`  

### Installation

Download the latest version of the software: https://github.com/snayfach/PhyloCNV/archive/v0.0.2.tar.gz 

Unpack the project: `tar -zxvf PhyloCNV-0.0.2.tar.gz`

Run setup.py. This will install any dependencies:  
`python setup.py install` or  
`sudo python setup.py install` to install as a superuser

Alternatively, you can manually install the software.
First, check that required python libraries are installed. You should be able to enter the following command in the python interpreter without getting an error:  
`>>> import Bio.SeqIO`  
`>>> import numpy`  
`>>> import pysam`  

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
  -d, --debug
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
  --species_profile     Estimate genome-cluster abundance
  --pangenome_build_db  Build bowtie2 database of pangenome centroids
  --pangenome_align     Align reads to genome-clusters
  --pangenome_cov       Compute coverage of pangenomes
  --snps_build_db       Build bowtie2 database of representative genomes
  --snps_align          Align reads to representative genomes
  --snps_call           Run samtools mpileup & estimate SNP frequencies

Species Abundance:
  --reads_gc READS_MS   # reads to use for estimating genome-cluster abundance
                        (5000000)

Species selection (choose one):
  --gc_topn GC_TOPN     Top N most abundant (None)
  --gc_cov GC_COV       Coverage threshold (None)
  --gc_rbun GC_RBUN     Relative abundance threshold (None)
  --gc_id GC_ID         Identifier of specific genome cluster or comma-
                        separated list of ids (None)

Pangenome module:
  --pangenome_align_speed {very-fast-local,fast-local,sensitive-local,very-sensitive-local}
                        alignment speed/sensitivity (very-sensitive-local)
  --pangenome_reads PANGENOME_READS
                        # reads for pangenome or genome alignment (use all)
  --pangenome_map_pid PANGENOME_MAP_PID
                        Minimum percent ID between read and reference (93.0)
  --pangenome_pid {90,92.5,95,97.5,99}
                        Reference gene cluster percent ID (97.5)
  --pangenome_aln_cov PANGENOME_ALN_COV
                        Minimum alignment coverage of read (0.70)

SNPs module:
  --snps_align_speed {very-fast,fast,sensitive,very-sensitive}
                        alignment speed/sensitivity (very-sensitive)
  --snps_reads SNPS_READS
                        # reads for pangenome or genome alignment (use all)
  --snps_mapq SNPS_MAPQ
                        Minimum map quality (20)
  --snps_baseq SNPS_BASEQ
                        Minimum base quality (20)  
```

### Output
Final outputs:  
* genome_clusters.abundance: abundance of 5,952 genome-clusters
* coverage: coverage of pan genes within each genome-cluster  
* snps: coverage, allele frequencies, and consensus alleles   

Intermediate files:  
* db: bowtie2 indexes of pangenomes and/or genomes   
* pangenome.bam: alignments against pangenome database 
* genomes.bam: alignments against genome database   
* genomes.vcf: vcf file generated from genomes.bam

### Example
Run PhyloCNV using the test FASTQ file and the top 5 most abundant species:
```
run_phylo_cnv.py \
-1 phylo_cnv/example/example.fastq.gz \
-D /path/to/genome_clusters \
-o phylo_cnv/example \
--all \
--gc_topn 5 \
--verbose
```

### Citation
If you use this tool, please cite:
Nayfach, S. and Pollard, KS. PhyloCNV: an integrated, high-resolution pipeline for quantifying strain-level variation from shotgun metagenomes (In Preparation).
