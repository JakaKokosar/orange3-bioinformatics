import io
import gzip
import csv

from requests.models import Response

from orangecontrib.bioinformatics.ncbi.taxonomy import common_taxids
from server_update.utils import RemoteResource


class HomologGenes(RemoteResource):

    def __init__(self):
        super().__init__()

        self.domain = 'homologene'
        self.file_name = 'homologene.tab'
        self.title = 'Homolog Genes'
        self.description = 'Homoloh genes from NCBI HomoloGene database.'
        self.tags = ['NCBI', 'genes', 'homologs', 'homolog genes', 'homologene']

        self.download_url = 'http://ftp.ncbi.nlm.nih.gov/pub/HomoloGene/build68/homologene.data'

    def load_gene_history(self):
        taxonomy_col = 0
        current_id_col = 1
        discontinued_id_col = 2

        response = self.fetch_url('http://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz')
        with gzip.open(io.BytesIO(response.content), 'rb') as file:
            # skip header
            file.readline()

            # store lines in memory
            gene_history = {}
            for line in file:
                info = tuple(line.decode().strip().split('\t'))
                gene_history[(info[taxonomy_col], info[discontinued_id_col])] = info[current_id_col]

        return gene_history

    def handle_response(self, response: Response) -> str:
        return response.content.decode('utf-8').strip()

    def prepare_data(self):
        response = self.fetch_url(self.download_url)
        homologs = self.handle_response(response)
        gene_history = self.load_gene_history()
        common_tax_ids = common_taxids()

        group_id = 0
        tax_id = 1
        entrez_id = 2

        with open('log.txt', 'a') as fp_log:

            homologs_out = [['group_id', 'tax_id', 'entrez_id']]
            for line in homologs.split('\n'):
                columns = line.split('\t')

                # Note: if '-' -> discontinued
                #       if 'some id' -> new entrez id
                #       if None -> valid entrez id
                gene_id_status = gene_history.get((columns[tax_id], columns[entrez_id]), None)

                if not columns[tax_id] in common_tax_ids:
                    continue
                elif gene_id_status == '-':
                    print(f'Gene {columns[entrez_id]} is discontinued', file=fp_log)
                    continue

                if gene_id_status is None:
                    gene_id = columns[entrez_id]
                else:
                    print(f'Gene {columns[entrez_id]} is now {gene_id_status}', file=fp_log)
                    gene_id = gene_id_status

                homologs_out.append([columns[group_id], columns[tax_id], gene_id])

        with open('homologene.tab', 'w') as fp:
            csv_writer = csv.writer(fp, delimiter='\t')
            csv_writer.writerows(homologs_out)


if __name__ == '__main__':
    homolog_db = HomologGenes()
    homolog_db.prepare_data()
    homolog_db.to_serverfile_format()
