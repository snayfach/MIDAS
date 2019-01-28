# MIDAS-IGGdb temporary documentation

The MIDAS-IGGdb tool under construction in this repo is based on the
original [MIDAS tool](https://github.com/snayfach/MIDAS) with modifications
that integrate all 23,790 species from [IGGdb](https://github.com/snayfach/IGGdb).

The 23,790 IGG species are derived by clustering a total of 206,581 genomes from
PATRIC, IMG, and the Human Gut MAG dataset to 95% identity.  In contrast, the original
MIDAS db version 1.2 represents only 30,000 genomes clustered to 5,926 species
at 99% identity.

A verbatim copy of [the original documentation](old-docs/README.md) is still
maintained in this repo while it is under construction.


# Species abundance estimation

To support the larger database, the original MIDAS `run_species` command
is being deprecated in favor of the new [IGGSearch tool](https://github.com/snayfach/IGGsearch).


# Building the MIDAS-IGG database

To accommodate the vast increase in source data and support its further enrichment
for SNP calling, we've created a MIDAS-IGG tool called `SMELTER` that supercedes
the legacy `build_midas_db` tool.

All inputs and outputs of SMELTER on the IGG data set are cached on the CZ Biohub
S3 under [s3://microbiome-bdimitrov/IGGdb](http://microbiome-bdimitrov.s3.amazonaws.com/IGGdb/README.TXT)
and open to public download.

SMELTER supports the following workflows.

## Building GMAP/GSNAP index for haplotype encrichment
```
cd /path/to/fast-scratch-space
/path/to/MIDAS-IGG/scripts/smelt.py build_gsnap_index_into outdir from /path/to/IGGdb/v1.0.0/metadata/species_info.tsv
mv outdir /path/to/IGGdb/v1.0.0/gsnap_repgenome_index
```

