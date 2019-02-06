#!/usr/bin/env python3
import random
import time
import json
from os.path import isfile, isdir, dirname, abspath, basename
from utilities import tsv_rows, parse_table, tsprint


class IGGdb:
    """Encapsulates the DB and provides an interface to look up information."""

    def __init__(self, iggdb_toc_species, quiet=False):
        try:
            assert basename(iggdb_toc_species) == "species_info.tsv"
        except Exception as e:
            e.help_text = f"Expected /path/to/species_info.tsv, was given '{iggdb_toc_species}' instead."
            raise
        try:
            self.iggdb_root = dirname(dirname(abspath(iggdb_toc_species)))
            iggdb_toc_genomes = f"{self.iggdb_root}/metadata/genome_info.tsv"
            assert isfile(iggdb_toc_genomes)
            assert isdir(f"{self.iggdb_root}/pangenomes")
            assert isdir(f"{self.iggdb_root}/repgenomes")
        except Exception as e:
            e.help_text = f"Unexpected MIDAS-IGGdb directory structure around {iggdb_toc_species}."
            raise
        self.species_info = list(parse_table(tsv_rows(iggdb_toc_species)))
        self.genome_info = list(parse_table(tsv_rows(iggdb_toc_genomes)))
        if not quiet:
            tsprint(f"Found {len(self.genome_info)} genomes in {iggdb_toc_genomes}.")
            tsprint(f"Found {len(self.species_info)} species in {iggdb_toc_species}, for example:")
            random.seed(time.time())
            random_index = random.randrange(0, len(self.species_info))
            tsprint(json.dumps(self.species_info[random_index], indent=4))
        self.species = {s['species_id']: s for s in self.species_info}
        self.genomes = {g['genome_id']: g for g in self.genome_info}
        for s in self.species_info:
            genome_id = s['representative_genome']
            g = self.genomes[genome_id]
            s['repgenome_with_origin'] = genome_id + "." + g['repository'].lower()
            s['repgenome_path'] = f"{self.iggdb_root}/repgenomes/{s['repgenome_with_origin']}.fna"
            s['pangenome_path'] = f"{self.iggdb_root}/pangenomes/{s['species_alt_id']}/centroids.fa"


    def get_species(self, species_id, default=None):
        return self.species.get(species_id, default)
