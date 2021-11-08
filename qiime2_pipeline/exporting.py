from .template import Processor, Settings


class ExportFeatureTable(Processor):

    feature_table_qza: str
    output_tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_table_qza: str, output_tsv: str):
        self.feature_table_qza = feature_table_qza
        self.output_tsv = output_tsv

        self.qza_to_biom()
        self.biom_to_tsv()

    def qza_to_biom(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.feature_table_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def biom_to_tsv(self):
        lines = [
            'biom convert --to-tsv',
            f'-i {self.workdir}/feature-table.biom',
            f'-o {self.output_tsv}'
        ]
        self.call(' \\\n  '.join(lines))


class ExportFeatureSequence(Processor):

    feature_sequence_qza: str
    output_fa: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            feature_sequence_qza: str,
            output_fa: str):

        self.feature_sequence_qza = feature_sequence_qza
        self.output_fa = output_fa

        self.qza_to_fa()
        self.move_fa()

    def qza_to_fa(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.feature_sequence_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def move_fa(self):
        self.call(f'mv {self.workdir}/dna-sequences.fasta {self.output_fa}')


class ExportAlignedSequence(Processor):

    aligned_sequence_qza: str
    output_fa: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            aligned_sequence_qza: str,
            output_fa: str):

        self.aligned_sequence_qza = aligned_sequence_qza
        self.output_fa = output_fa

        self.qza_to_fa()
        self.move_fa()

    def qza_to_fa(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.aligned_sequence_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def move_fa(self):
        self.call(f'mv {self.workdir}/aligned-dna-sequences.fasta {self.output_fa}')


class ExportTree(Processor):

    tree_qza: str
    output_nwk: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            tree_qza: str,
            output_nwk: str):

        self.tree_qza = tree_qza
        self.output_nwk = output_nwk

        self.qza_to_nwk()
        self.move_nwk()

    def qza_to_nwk(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.tree_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def move_nwk(self):
        self.call(f'mv {self.workdir}/tree.nwk {self.output_nwk}')


class ExportTaxonomy(Processor):

    taxonomy_qza: str
    output_tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, taxonomy_qza: str, output_tsv: str):

        self.taxonomy_qza = taxonomy_qza
        self.output_tsv = output_tsv

        self.qza_to_tsv()
        self.move_tsv()

    def qza_to_tsv(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.taxonomy_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def move_tsv(self):
        self.call(f'mv {self.workdir}/taxonomy.tsv {self.output_tsv}')
