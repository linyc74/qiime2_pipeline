from os.path import basename
from .tools import edit_fpath
from .template import Processor, Settings


class ExportFeatureTable(Processor):

    feature_table_qza: str
    tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_table_qza: str) -> str:
        self.feature_table_qza = feature_table_qza

        self.qza_to_biom()
        self.biom_to_tsv()

        return self.tsv

    def qza_to_biom(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.feature_table_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def biom_to_tsv(self):
        self.tsv = edit_fpath(
            fpath=self.feature_table_qza,
            old_suffix='.qza',
            new_suffix='.tsv',
            dstdir=self.workdir
        )
        lines = [
            'biom convert --to-tsv',
            f'-i {self.workdir}/feature-table.biom',
            f'-o {self.tsv}'
        ]
        self.call(' \\\n  '.join(lines))


class ExportFeatureSequence(Processor):

    feature_sequence_qza: str
    output_fa: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_sequence_qza: str) -> str:

        self.feature_sequence_qza = feature_sequence_qza

        self.qza_to_fa()
        self.move_fa()

        return self.output_fa

    def qza_to_fa(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.feature_sequence_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def move_fa(self):
        self.output_fa = edit_fpath(
            fpath=self.feature_sequence_qza,
            old_suffix='.qza',
            new_suffix='.fa',
            dstdir=self.workdir
        )
        self.call(f'mv {self.workdir}/dna-sequences.fasta {self.output_fa}')


class ExportTaxonomy(Processor):

    taxonomy_qza: str
    tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, taxonomy_qza: str) -> str:
        self.taxonomy_qza = taxonomy_qza

        self.qza_to_tsv()
        self.move_tsv()

        return self.tsv

    def qza_to_tsv(self):
        self.tsv = edit_fpath(
            fpath=self.taxonomy_qza,
            old_suffix='.qza',
            new_suffix='.tsv',
            dstdir=self.workdir
        )
        lines = [
            'qiime tools export',
            f'--input-path {self.taxonomy_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def move_tsv(self):
        self.call(f'mv {self.workdir}/taxonomy.tsv {self.tsv}')
