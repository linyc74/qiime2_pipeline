from os.path import abspath
from .tools import edit_fpath
from .template import Processor, Settings


class Export(Processor):

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def mv(self, src: str, dst: str):
        if abspath(src) != abspath(dst):
            self.call(f'mv {src} {dst}')


class ExportFeatureTable(Export):

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
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools export',
            f'--input-path {self.feature_table_qza}',
            f'--output-path {self.workdir}',
        ])
        self.call(cmd)

    def biom_to_tsv(self):
        self.tsv = edit_fpath(
            fpath=self.feature_table_qza,
            old_suffix='.qza',
            new_suffix='.tsv',
            dstdir=self.workdir
        )
        cmd = self.CMD_LINEBREAK.join([
            'biom convert --to-tsv',
            f'-i {self.workdir}/feature-table.biom',
            f'-o {self.tsv}'
        ])
        self.call(cmd)


class ExportFeatureSequence(Export):

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
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools export',
            f'--input-path {self.feature_sequence_qza}',
            f'--output-path {self.workdir}',
        ])
        self.call(cmd)

    def move_fa(self):
        self.output_fa = edit_fpath(
            fpath=self.feature_sequence_qza,
            old_suffix='.qza',
            new_suffix='.fa',
            dstdir=self.workdir
        )
        self.mv(f'{self.workdir}/dna-sequences.fasta', self.output_fa)


class ExportTaxonomy(Export):

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
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools export',
            f'--input-path {self.taxonomy_qza}',
            f'--output-path {self.workdir}',
        ])
        self.call(cmd)

    def move_tsv(self):
        self.tsv = edit_fpath(
            fpath=self.taxonomy_qza,
            old_suffix='.qza',
            new_suffix='.tsv',
            dstdir=self.workdir
        )
        self.mv(f'{self.workdir}/taxonomy.tsv', self.tsv)


class ExportAlignedSequence(Export):

    aligned_sequence_qza: str
    output_fa: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, aligned_sequence_qza: str) -> str:
        self.aligned_sequence_qza = aligned_sequence_qza

        self.qza_to_fa()
        self.move_fa()

        return self.output_fa

    def qza_to_fa(self):
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools export',
            f'--input-path {self.aligned_sequence_qza}',
            f'--output-path {self.workdir}',
        ])
        self.call(cmd)

    def move_fa(self):
        self.output_fa = edit_fpath(
            fpath=self.aligned_sequence_qza,
            old_suffix='.qza',
            new_suffix='.fa',
            dstdir=self.workdir
        )
        self.mv(f'{self.workdir}/aligned-dna-sequences.fasta', self.output_fa)


class ExportTree(Export):

    tree_qza: str
    nwk: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, tree_qza: str) -> str:
        self.tree_qza = tree_qza

        self.qza_to_nwk()
        self.move_nwk()

        return self.nwk

    def qza_to_nwk(self):
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools export',
            f'--input-path {self.tree_qza}',
            f'--output-path {self.workdir}',
        ])
        self.call(cmd)

    def move_nwk(self):
        self.nwk = edit_fpath(
            fpath=self.tree_qza,
            old_suffix='.qza',
            new_suffix='.nwk',
            dstdir=self.workdir
        )
        self.mv(f'{self.workdir}/tree.nwk', self.nwk)


class ExportBetaDiversity(Export):

    distance_matrix_qza: str
    tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, distance_matrix_qza: str) -> str:
        self.distance_matrix_qza = distance_matrix_qza
        self.qza_to_tsv()
        self.move_tsv()
        return self.tsv

    def qza_to_tsv(self):
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools export',
            f'--input-path {self.distance_matrix_qza}',
            f'--output-path {self.workdir}',
        ])
        self.call(cmd)

    def move_tsv(self):
        self.tsv = edit_fpath(
            fpath=self.distance_matrix_qza,
            old_suffix='.qza',
            new_suffix='.tsv',
            dstdir=self.workdir
        )
        self.mv(f'{self.workdir}/distance-matrix.tsv', self.tsv)
