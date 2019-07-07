import gzip

from io import BytesIO
from collections import defaultdict
from requests.models import Response

from Orange.data import Table, Domain, StringVariable
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher, Gene
from server_update.utils import RemoteResource


class PanglaoDB(RemoteResource):

    def __init__(self):
        super().__init__()

        self.domain = 'marker_genes'
        self.file_name = 'panglao_gene_markers.tab'
        self.title = 'Marker Genes: PanglaoDB'
        self.description = 'Marker genes from Panglao database.'
        self.tags = ['panglao', 'genes', 'markers', 'marker genes']

        self.reference, self.reference_url = 'PanglaoDB', 'https://panglaodb.se/'
        self.download_url = 'https://panglaodb.se/markers/PanglaoDB_markers_15_May_2019.tsv.gz'

    def handle_response(self, response: Response):
        with gzip.GzipFile(fileobj=BytesIO(response.content), mode='r') as f:
            return f.read().decode('utf-8').strip()

    def prepare_data(self) -> Table:
        response = self.fetch_url(self.download_url)
        content = self.handle_response(response)

        species = 0
        gene_symbol = 1
        cell_type = 2
        genes_by_organism = defaultdict(list)
        organism_mapper = {'Mm': 'Mouse', 'Hs': 'Human'}

        def _gene_function_table(desc_col: StringVariable, gm_results: GeneMatcher):
            _domain = Domain([], metas=[desc_col])
            _data = [[str(gene.description) if gene.description else ''] for gene in gm_results.genes]
            return Table(_domain, _data)

        for line in content.split('\n'):
            columns = line.split('\t')

            for org in columns[species].split(' '):
                if org in organism_mapper.keys():
                    gene_entry = [organism_mapper[org],
                                  columns[gene_symbol],
                                  columns[cell_type],
                                  self.reference,
                                  self.reference_url]
                    genes_by_organism[organism_mapper[org]].append(gene_entry)

        domain = Domain([], metas=[StringVariable('Organism'),
                                   StringVariable('Name'),
                                   StringVariable('Cell Type'),
                                   StringVariable('Reference'),
                                   StringVariable('URL')])

        entrez_id_column = StringVariable('Entrez ID')
        description_column = StringVariable('Function')

        # construct data table for mouse
        gm_mouse = GeneMatcher('10090', case_insensitive=True)
        mouse_table = Table(domain, genes_by_organism['Mouse'])
        mouse_table = gm_mouse.match_table_column(mouse_table, 'Name', entrez_id_column)
        mouse_table = Table.concatenate([mouse_table, _gene_function_table(description_column, gm_mouse)])

        # construct data table for human
        gm_human = GeneMatcher('9606', case_insensitive=True)
        human_table = Table(domain, genes_by_organism['Human'])
        human_table = gm_human.match_table_column(human_table, 'Name', entrez_id_column)
        human_table = Table.concatenate([human_table, _gene_function_table(description_column, gm_human)])

        # return combined tables
        return Table.concatenate([mouse_table, human_table], axis=0)


class CellMarkerDB(RemoteResource):

    def __init__(self):
        super().__init__()

        self.domain = 'marker_genes'
        self.file_name = 'cellMarker_gene_markers.tab'
        self.title = 'Marker Genes: CellMarker'
        self.description = 'Marker genes from CellMarker database.'
        self.tags = ['CellMarker', 'genes', 'markers', 'marker genes']

        self.reference, self.reference_url = '{}', 'https://www.ncbi.nlm.nih.gov/pubmed/?term={}'
        self.download_url = 'http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/all_cell_markers.txt'

    def handle_response(self, response: Response):
        return response.content.decode('utf-8').strip()

    def prepare_data(self) -> Table:
        response = self.fetch_url(self.download_url)
        content = self.handle_response(response)

        species_col = 0
        cell_name_col = 5
        gene_symbol_col = 8
        entrez_id_col = 9
        pubmed_id_col = 13
        organisms = ['Human', 'Mouse']
        data = list()

        domain = Domain([], metas=[StringVariable('Organism'),
                                   StringVariable('Name'),
                                   StringVariable('Cell Type'),
                                   StringVariable('Reference'),
                                   StringVariable('URL'),
                                   StringVariable('Entrez ID')])

        unique_rows = set()

        for line in content.split('\n'):
            columns = line.split('\t')
            organism = columns[species_col]

            if organism in organisms:
                symbols = columns[gene_symbol_col].replace('[', '').replace(']', '').split(', ')
                entrez_ids = columns[entrez_id_col].replace('[', '').replace(']', '').split(', ')

                for symbol, entrez_id in zip(symbols, entrez_ids):
                    try:
                        int(entrez_id)
                    except ValueError:
                        continue

                    if (columns[cell_name_col], entrez_id) in unique_rows:
                        continue

                    ref = self.reference.format(columns[pubmed_id_col])

                    try:
                        int(ref)
                        ref_link = self.reference_url.format(columns[pubmed_id_col])
                    except ValueError:
                        # its not pubmed_id
                        ref_link = '?'

                    unique_rows.add((columns[cell_name_col], entrez_id))
                    gene_entry = [organism,
                                  symbol,
                                  columns[cell_name_col],
                                  ref,
                                  ref_link,
                                  entrez_id]

                    data.append(gene_entry)

        table = Table(domain, data)
        genes = [Gene(name) for name in table.get_column_view('Entrez ID')[0]]
        for gene in genes:
            gene.ncbi_id = gene.input_name
            gene.load_ncbi_info()

        description_column = StringVariable('Function')
        domain = Domain([], metas=table.domain.metas + (description_column, ))
        table = table.transform(domain)
        table[:, description_column] = [[gene.description] for gene in genes]
        return table


if __name__ == '__main__':
    pang_db = PanglaoDB()
    data_table = pang_db.prepare_data()
    data_table.save(pang_db.file_name)
    pang_db.to_serverfile_format()

    cell_marker_db = CellMarkerDB()
    data_table = cell_marker_db.prepare_data()
    data_table.save(cell_marker_db.file_name)
    cell_marker_db.to_serverfile_format()

    # from serverfiles import ServerFiles
    #
    # with open('__INFO__', 'wt') as f:
    #     # we must initialize ServerFiles object again because old one has __INFO__ cached
    #     json.dump(list(ServerFiles(server='http://download.biolab.si/datasets/bioinformatics/').allinfo().items()), f)

    # serverfiles = ServerFiles(server='http://download.biolab.si/datasets/bioinformatics/')
    # from orangecontrib.bioinformatics.utils import serverfiles
    #
    # serverfiles.update('marker_genes', 'panglao_gene_markers.tab')
    # file_path = serverfiles.localpath_download('marker_genes', 'panglao_gene_markers.tab')
    # print(file_path)
    #
