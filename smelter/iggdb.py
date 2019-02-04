#!/usr/bin/env python3
from smelter.utilities import backtick, tsv_rows, parse_table


def repgenome_for(filename, #pylint: disable=dangerous-default-value
                  ORIGINS=["patric", "hgm", "img"],
                  EXTENSIONS=["fna"]):
    """Deconstructs a repgenome filename into constituent elements, which can be used to look up information in IGGdb species_info"""
    rg = {}
    rg['id'], rg['origin'], rg['extension'] = filename.rsplit('.', 2)
    assert rg['extension'] in EXTENSIONS
    assert rg['origin'] in ORIGINS
    rg['id.origin'] = rg['id'] + '.' + rg['origin']
    rg['id.origin.extension'] = rg['id.origin'] + '.' + rg['extension']
    return rg


class IGGdb:
    """Encapsulates the DB and provides an interface to look up information."""

    def __init__(self, iggdb_root):
        self.iggdb_root = iggdb_root
        self.species = list(parse_table(tsv_rows(f"{iggdb_root}/metadata/species_info.tsv")))
        # Glob the repgenome directory and populate species info with corresponding repgenome paths.
        repgenome_references = {}
        for filename in backtick(f"ls {iggdb_root}/repgenomes").strip().split("\n"):
            rg = repgenome_for(filename)
            repgenome_references[rg['id']] = rg
        for s in self.species:
            s_rg = repgenome_references[s['representative_genome']]
            s['representative_genome_with_origin'] = s_rg['id.origin']
            s['representative_genome_path'] = f"{iggdb_root}/repgenomes/{s_rg['id.origin.extension']}"
        self.species_by_id = {}
        for s in self.species:
            self.species_by_id[s['species_id']] = s

    def get_species(self, species_id, default=None):
        return self.species_by_id.get(species_id, default)
