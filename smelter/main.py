#!/usr/bin/env python3
#
# MIDAS IGG SMELTER
#
# A tool to repackage and enrich information provided by Stephen Nayfach
# via http://github.com/snayfach/IGGdb in order to support genotyping and
# SNP calling over the IGGdb dataset.
#
# Created on January 28, 2019 by Boris Dimitrov, bdimitrov@chanzuckerberg.com.
# Distributed freely under the MIT License, subject to restrictions set forth
# in the licensing terms of MIDAS and IGG tools by Stephen Nayfach.  Please
# see http://github.com/snayfach/MIDAS and http://github.com/snayfach/IGGsearch
# for appropriate citation requirements.

import sys
assert sys.version_info >= (3, 6), "Please run this script with Python version >= 3.6."


import os
import traceback
import json
import time
from utilities import tsprint, backtick, makedirs, ProgressTracker
from iggdb import IGGdb


def smelt(argv):
    cwd = backtick('pwd')
    my_command = f"cd {cwd}; " + ' '.join(argv)
    tsprint(my_command)
    _, subcmd, outdir, iggdb_toc = argv
    subcmd = subcmd.replace("-", "_").lower()
    SUBCOMMANDS = {
        f"collate_{gdim}": gdim
        for gdim in ["pangenomes", "repgenomes"]
    }
    gdim = SUBCOMMANDS.get(subcmd)
    try:
        assert gdim in ["pangenomes", "repgenomes"]
    except Exception as e:
        e.help_text = f"Try a supported subcommand instead of {subcmd}."
        raise
    makedirs(outdir, exist_ok=False)
    iggdb = IGGdb(iggdb_toc)
    tsprint(f"Now collating fasta for gsnap {gdim} index construction.")
    MAX_FAILURES = 100
    count_successes = 0
    ticker = ProgressTracker(target=len(iggdb.species_info))
    failures = []
    for s in iggdb.species_info:
        try:
            s_species_alt_id = s['species_alt_id']
            s_tempfile = f"{outdir}/temp_{gdim}_{s_species_alt_id}.fa"
            # Note how the header tags we emit below for pangenomes and repgenomes are consistent;
            # this should enable easy reconciliation of gsnap alignments against the
            # two separate indexes.
            if gdim == "pangenomes":
                #
                # The header tag we wish to emit would be
                #
                #    >4547837|1657.8.patric|rest_of_original_header_from_pangenome_file
                #
                # where
                #
                #    species_alt_id = 4547837                  # from species_alt_id column in table
                #    repgenome_with_origin = 1657.8.patric     # from original header in pangenome file
                #
                # As the original header in the pangenome file already begins
                # with s_rg_id_origin, we just need to prepend species_alt_id.
                #
                s_header_xform = f"sed 's=^>=>{s_species_alt_id}|=' {s['pangenome_path']} > {s_tempfile} && cat {s_tempfile} >> {outdir}/temp_{gdim}.fa && rm {s_tempfile} && echo SUCCEEDED || echo FAILED"
            else:
                assert gdim == "repgenomes"
                #
                # The header tag we wish to emit would be
                #
                #    >4547837|1657.8.patric|entire_original_header_from_repgenome_file
                #
                # where
                #
                #    species_alt_id = 4547837                # species_alt_id column in species_info
                #    repgenome_with_origin = 1657.8.patric   # from file listing in repgenomes dir
                #
                s_repgenome_with_origin = s['repgenome_with_origin']
                s_repgenome_path = s['repgenome_path']
                s_header_xform = f"sed 's=^>=>{s_species_alt_id}|{s_repgenome_with_origin}|=' {s_repgenome_path} > {s_tempfile} && cat {s_tempfile} >> {outdir}/temp_{gdim}.fa && rm {s_tempfile} && echo SUCCEEDED || echo FAILED"
            status = backtick(s_header_xform)
            assert status == "SUCCEEDED"
            count_successes += 1
        except Exception as e:
            failures.append(s)
            if len(failures) == MAX_FAILURES:
                count_examined = len(failures) + count_successes
                e.help_text = f"Giving up after {MAX_FAILURES} failures in first {count_examined} species.  See temp files for more info."
                raise
        finally:
            ticker.advance(1)
    failed_species_alt_ids = [s['species_alt_id'] for s in failures]
    if not failures:
        tsprint(f"All {len(iggdb.species_info)} species were processed successfully.")
    else:
        tsprint(f"Collation of {len(failures)} species failed.  Those are missing from the final {gdim}.fa")
    # Create output file only on success.
    # Dump stats in json.
    collation_status = {
        "comment": f"Collation into {gdim}.fa succeeded on {time.asctime()} with command '{my_command}'.",
        "successfully_collated_species_count": count_successes,
        "failed_species_count": len(failures),
        "total_species_count": len(iggdb.species_info),
        "failed_species_alt_ids": failed_species_alt_ids,
        "elapsed_time": ticker.t_elapsed
    }
    collation_status_str = json.dumps(collation_status, indent=4)
    with open(f"{outdir}/{gdim}_collation_status.json", "w") as pcs:
        chars_written = pcs.write(collation_status_str)
        assert chars_written == len(collation_status_str)
        tsprint(collation_status_str)
    os.rename(f"{outdir}/temp_{gdim}.fa", f"{outdir}/{gdim}.fa")

def main():
    try:
        smelt(sys.argv)
    except Exception as e:
        tsprint(traceback.format_exc())
        tsprint("*** USAGE:  See https://github.com/czbiohub/MIDAS-IGGdb/blob/master/README.md#smelter ***\n")
        if hasattr(e, 'help_text'):
            tsprint(f"*** {e.help_text} ***") # pylint: disable=no-member


if __name__ == "__main__":
    main()
