# Temporary documentation

The MIDAS-IGGdb tool under construction in this repo is based on the
original [MIDAS tool](https://github.com/snayfach/MIDAS) with modifications
that integrate all 23,790 species from [IGGdb](https://github.com/snayfach/IGGdb).
The original MIDAS db version 1.2 spans only 5,926 species.

To support the larger database, the original MIDAS `run_species` command
is being deprecated in favor of the new [IGGSearch tool](https://github.com/snayfach/IGGsearch).
Details on how this affects various workflows will follow shortly.

A verbatim copy of [the original documentation](old-docs/README.md) is still
maintained in this repo while it is under construction.

In total, this IGGdb update represents 206,581 genomes fromm PATRIC, IMG, and
the Human Gut MAG dataset clustered to 23,790 species at 95% identity.  All
data required to construct the new MIDAS-IGGdb database is currently cached
on the CZ Biohub S3 under [s3://microbiome-bdimitrov/IGGdb](http://microbiome-bdimitrov.s3.amazonaws.com/IGGdb/README.TXT)
and completely open to public download.
