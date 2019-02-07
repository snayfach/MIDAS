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


# Species abundance estimation

To support the larger database, the original MIDAS `run_species` command
is being deprecated in favor of the new [IGGSearch tool](https://github.com/snayfach/IGGsearch).

# Variant calling

Here is how to run the original MIDAS bowtie2 based workflow for counting SNPs, on the larger IGGdb database.

1) Download `s3://microbiome-bdimitrov/IGGdb/v1.0.0/{metadata, iggsearch, repgenomes}` and clone `github.com/czbiohub/MIDAS-IGGdb` right next to each other

2) For any `sample_dir` of your choice, `mkdir -p sample_dir/results/iggsearch`

3) run `iggsearch` as per its own usage instructions, with output dir `sample_dir/results/iggsearch`

4) `/path/to/MIDAS-IGGdb/run.sh snps sample_dir/results -1 sample_dir/*.1.fastq -2 sample_dir/*.2.fastq -t 24` assuming a machine with 24 physical CPU cores.

Results will appear in `sample_dir/results/snps`.


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
