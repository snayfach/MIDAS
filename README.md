# Temporary documentation

The MIDAS-IGGdb tool under construction in this repo is based on the
original [MIDAS tool](https://github.com/snayfach/MIDAS) with modifications
that integrate all 23,790 species from [IGGdb](https://github.com/snayfach/IGGdb).
The original MIDAS db version 1.2 spans only 5,926 species.

To support the larger database, the original MIDAS `run_species` command
is being deprecated in favor of the new [IGGSearch tool](https://github.com/snayfach/IGGsearch).
A new tool called [Smelter](#smelter) helps repackage and refine IGG data
to support the various MIDAS-IGG workflows.

A verbatim copy of [the original documentation](old-docs/README.md) is still
maintained in this repo while it is under construction.

In total, this IGGdb update represents 206,581 genomes fromm PATRIC, IMG, and
the Human Gut MAG dataset clustered to 23,790 species at 95% identity.  All
data required to construct the new MIDAS-IGGdb database is currently cached
on the CZ Biohub S3 under [s3://microbiome-bdimitrov/IGGdb](http://microbiome-bdimitrov.s3.amazonaws.com/IGGdb/README.TXT)
and completely open to public download.

# SMELTER

To accommodate the vast increase in source data and support further enrichment
for SNP calling, we've built a new tool called `SMELTER` that supercedes
the legacy `build_midas_db` tool.

## Construct GMAP/GSNAP index for haplotype encrichment

```
cd /path/to/fast-scratch-space
/path/to/MIDAS-IGG/scripts/smelt.py build_gsnap_index candiate_gsnap_index /path/to/IGGdb/v1.0.0/metadata/species_info.tsv
mv candidate_gsnap_index /path/to/IGGdb/v1.0.0/gsnap_index
```
