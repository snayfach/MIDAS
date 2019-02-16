# MIDAS-IGGdb temporary documentation

The MIDAS-IGGdb tool under construction in this repo is based on the
original [MIDAS tool](https://github.com/snayfach/MIDAS) with modifications
that integrate all 23,790 species from [IGGdb](https://github.com/snayfach/IGGdb),
which are derived by clustering a total of 206,581 genomes from
PATRIC, IMG, and the Human Gut MAG dataset to 95% identity.

In contrast, the original
MIDAS db version 1.2 represents only 30,000 genomes clustered to 5,926 species
at 99% identity.

A verbatim copy of [the original documentation](old-docs/README.md) is still
maintained in this repo while it is under construction.

# Setup

Download `s3://microbiome-bdimitrov/IGGdb/v1.0.0/{metadata, iggsearch, repgenomes}` and clone `github.com/czbiohub/MIDAS-IGGdb` right next to each other.  This lets the `run.sh` script find everything automatically.

# Species abundance estimation

To support the larger database, the original MIDAS `run_species` command
is being deprecated in favor of the new [IGGSearch tool](https://github.com/snayfach/IGGsearch).

# Variant calling - single sample

1) For any `sample_dir` of your choice, run `iggsearch` with output dir `sample_dir/results/iggsearch`

2) `/path/to/MIDAS-IGGdb/run.sh snps sample_dir/results -1 sample_dir/*.1.fastq -2 sample_dir/*.2.fastq -t 24`

Results will appear in `sample_dir/results/snps`.

# Variant calling - pooled

After running the single-sample variant calling workflow above for several samples, results can be merged and annotated as follows.

```
OUTDIR=sample_dir_xyz/..
cd $OUTDIR
ls sample_dir_*/results > samples.txt   # list of all samples to run on
/path/to/MIDAS-IGGdb/run.sh merge snps . -i samples.txt -t file --threads 24
```
Results will appear in `${OUTDIR}/OTU-*`  (one folder per species present).

# *NEW* Variant calling - batch of samples - fixed subset of species

This workflow is new in MIDAS-IGGdb.  It builds a bowtie DB for a fixed set of species, then aligns a set of samples against that DB, and counts SNPS in each sample.  Genomic positions with zero mapped reads are omitted from the output.  Tested with 975 species and 40 samples.  Could easily do more.

1. Build a bowtie DB using the new `--all_species_in_db` and `--dbtoc` options to specify your desired subset of species.
```
/path/to/MIDAS-IGGdb/run.sh snps sample_dir/results -t 48 --build-db --all_species_in_db --dbtoc /path/to/IGGdb/v1.0.0/gtpro_subset_975/species_info.tsv
```

2. Save that DB to a central location.
```
mv sample_dir/results/temp /path/to/IGGdb/v1.0.0/gtpro_subset_975_bowtie2_repgenomes
```

3. Perform alignment on all samples.
```
for f in <....list of sample dirs....>; do
    test -e $f/results/snps/temp/genomes.bam && echo "Alignment already completed for $f" && continue
    rm -f $f.SUCCESS && touch $f.FAILURE
    /path/to/MIDAS-IGGdb/run.sh snps $f/results -1 $f/sample_reads.fastq -t 48 --align --all_species_in_db --dbtoc /path/to/IGGdb/v1.0.0/gtpro_subset_975/species_info.tsv --bowtie-db /path/to/IGGdb/v1.0.0/gtpro_subset_975_bowtie_repgenomes && rm $f.FAILURE && touch $f.SUCCESS
done
```

4. Perform pileup on all samples.
```
for f in <....list of sample dirs....>; do
    test -e $f/results/snps/summary.txt && echo "Pileup already completed for $f" && continue
    rm -f $f.SUCCESS && touch $f.FAILURE
    /path/to/MIDAS-IGGdb/run.sh snps $f/results -1 $f/sample_reads.fastq -t 48 --pileup --sparse --all_species_in_db --dbtoc /path/to/IGGdb/v1.0.0/gtpro_subset_975/species_info.tsv && rm $f.FAILURE && touch $f.SUCCESS
done
```
When using a remote machine, it's best to save the above snippets into script files and run under nohup or screen.

Each stage will emit SUCCESS or FAILURE files as appropriate for every sample.  The pileup for each genome in the specified subset will be in sample_dir/output (for each sample_dir).

# Building the MIDAS-IGGdb database

To accommodate the vast increase in source data and support its further enrichment
for SNP calling, we've created a MIDAS-IGG tool called `SMELTER` that supercedes
the legacy `build_midas_db` tool.

All inputs and outputs of SMELTER on the IGG data set are cached on the CZ Biohub
S3 under [s3://microbiome-bdimitrov/IGGdb](http://microbiome-bdimitrov.s3.amazonaws.com/IGGdb/README.TXT)
and open to public download.

SMELTER supports the following workflows.

## Construct GSNAP DB for all IGGdb genomes
```
/path/to/smelter/main.py collate_repgenomes /fast-scratch-space/new_work_dir /path/to/IGGdb/v1.0.0/metadata/species_info.tsv
cd /fast-scratch-space/new_work_dir
nohup time gmap_build -D . -d repgenomes_nt16 -k 16 repgenomes.fa
rmdir repgenomes_nt16.maps
mv /fast-scratch-space/new_work_dir /path/to/IGGdb/v1.0.0/gsnap_repgenomes
```
Then repeat for pangenomes instead of repgenomes.

This also works on a subset of the entire species_info.tsv file, e.g. `/path/to/IGGdb/v1.0.0/my_subset_metadata/species_info.tsv`
