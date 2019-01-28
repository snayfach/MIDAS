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
from utilities import tsprint, backtick, makedirs


def smelt(argv):
    tsprint(f"cd {backtick('pwd')}")
    tsprint(' '.join(argv))
    _, subcmd, outdir, iggtoc = argv
    assert subcmd.replace("-", "_").lower() == "build_gsnap_index"
    assert os.path.basename(iggtoc) == "species_info.tsv"
    makedirs(outdir, exist_ok=False)


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
