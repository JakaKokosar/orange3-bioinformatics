import sqlite3
import gzip
import json
import io
import requests_cache

from requests import Response
from typing import Set, Tuple, List
from collections import defaultdict
from sqlite_utils import Database

from orangecontrib.bioinformatics.ncbi.gene import parse_sources, parse_synonyms
from orangecontrib.bioinformatics.ncbi.taxonomy import Taxonomy, common_taxids, common_taxid_to_name
from server_update.utils import RemoteResource


requests_cache.install_cache('cache_gene_info')
taxonomy_db = Taxonomy()
relevant_taxonomy_ids = set()
gene_info = defaultdict(set)

# columns indexes
# ftp://ftp.ncbi.nlm.nih.gov/gene/README under "gene_info" section
tax_id, gene_id, symbol, synonyms, db_refs, description = 0, 1, 2, 4, 5, 8
locus_tag, chromosome, map_location, type_of_gene, modification_date = 3, 6, 7, 9, 14
symbol_from_nomenclature_authority, full_name_from_nomenclature_authority = 10, 11
nomenclature_status, other_designations = 12, 13

init_table = """
CREATE TABLE "gene_info" (
    species TEXT NOT NULL,
    tax_id TEXT NOT NULL,
    gene_id TEXT NOT NULL UNIQUE,
    symbol TEXT NOT NULL,
    synonyms TEXT,
    db_refs TEXT,
    description TEXT,
    locus_tag TEXT,
    chromosome TEXT,
    map_location TEXT,
    type_of_gene TEXT,
    symbol_from_nomenclature_authority TEXT,
    full_name_from_nomenclature_authority TEXT,
    nomenclature_status TEXT,
    other_designations TEXT,
    modification_date TEXT);
"""


class GeneInfo(RemoteResource):

    def __init__(self, taxonomy_id):
        super().__init__()

        self.taxonomy_id = taxonomy_id
        self.species_name = taxonomy_db.get_entry(taxonomy_id).name

        self.domain: str = 'gene_info'
        self.title: str = f'Gene Info: {self.species_name}'
        self.file_name: str = f'{taxonomy_id}.sqlite'
        self.description: str = f'Gene Info for {self.species_name} and all its strains.'
        self.tags: List[str] = ['NCBI', 'genes', 'info', 'gene info', self.species_name, taxonomy_id]

        self.download_url: str = 'http://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz'

    def handle_response(self, response: Response):
        global gene_info, relevant_taxonomy_ids

        if len(gene_info) != 0:
            return

        print('Parsing gene_info.gz ...')
        with gzip.open(io.BytesIO(response.content), 'rb') as info_file:
            # skip header
            info_file.readline()

            # store lines in memory
            gene_info = defaultdict(set)
            for line in info_file:
                info = tuple(line.decode().strip().split('\t'))

                if info[tax_id] in relevant_taxonomy_ids:
                    gene_info[info[tax_id]].add(info)

    def prepare_data(self):
        print(f'Creating database for {taxonomy_db.get_entry(self.taxonomy_id).name} ...')

        response = self.fetch_url(self.download_url)

        self.handle_response(response)

        species_genes: Set[Tuple] = set()
        for _tax in [self.taxonomy_id] + taxonomy_db.get_all_strains(self.taxonomy_id):
            species_genes.update(gene_info.get(_tax, []))

        con = sqlite3.connect(self.file_name, timeout=15)
        db = Database(con)

        cursor = con.cursor()

        cursor.execute("DROP TABLE IF EXISTS 'gene_info'")
        cursor.execute(init_table)

        cursor.executemany(
            "INSERT INTO gene_info VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                (
                    common_taxid_to_name(self.taxonomy_id),
                    gene[tax_id],
                    gene[gene_id],
                    gene[symbol],
                    json.dumps({'values': parse_synonyms(gene[synonyms])}),
                    json.dumps(parse_sources(gene[db_refs])),
                    gene[description] if gene[description] != '-' else None,
                    gene[locus_tag] if gene[locus_tag] != '-' else None,
                    gene[chromosome] if gene[chromosome] != '-' else None,
                    gene[map_location] if gene[map_location] != '-' else None,
                    gene[type_of_gene] if gene[type_of_gene] != '-' else None,
                    gene[symbol_from_nomenclature_authority] if gene[symbol_from_nomenclature_authority] != '-' else None,
                    gene[full_name_from_nomenclature_authority] if gene[full_name_from_nomenclature_authority] != '-' else None,
                    gene[nomenclature_status] if gene[nomenclature_status] != '-' else None,
                    gene[other_designations] if gene[other_designations] != '-' else None,
                    gene[modification_date],
                )
                for gene in species_genes
            ),
        )

        con.commit()
        db['gene_info'].enable_fts(['gene_id', 'symbol', 'synonyms', 'db_refs',
                                    'description', 'locus_tag', 'symbol_from_nomenclature_authority'])
        db['gene_info'].optimize()

        con.close()
        print()


if __name__ == '__main__':
    import time

    print('Loading taxonomies ...')
    for tax in common_taxids():
        relevant_taxonomy_ids.add(tax)
        relevant_taxonomy_ids.update(taxonomy_db.get_all_strains(tax))

    st = time.time()

    for species in common_taxids():
        gi = GeneInfo(species)
        gi.prepare_data()
        gi.to_serverfile_format()

    et = time.time()
    print(f'This took {et - st} seconds! ')

