from os.path import abspath
from .tools import edit_fpath
from .template import Processor


class Export(Processor):

    def qiime_tools_export(self, input_path: str, output_path: str):
        log = f'{self.outdir}/qiime-tools-export.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools export',
            f'--input-path {input_path}',
            f'--output-path {output_path}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)

    def mv(self, src: str, dst: str):
        if abspath(src) != abspath(dst):
            self.call(f'mv "{src}" "{dst}"')


class ExportFeatureTable(Export):

    feature_table_qza: str
    tsv: str

    def main(self, feature_table_qza: str) -> str:
        self.feature_table_qza = feature_table_qza

        self.qza_to_biom()
        self.biom_to_tsv()
        self.remove_first_line()  # remove the first line of the tsv file: "# Constructed from biom file"

        return self.tsv

    def qza_to_biom(self):
        self.qiime_tools_export(
            input_path=self.feature_table_qza,
            output_path=self.workdir)

    def biom_to_tsv(self):
        self.tsv = edit_fpath(
            fpath=self.feature_table_qza,
            old_suffix='.qza',
            new_suffix='.tsv',
            dstdir=self.workdir
        )
        log = f'{self.outdir}/biom-convert.log'
        cmd = self.CMD_LINEBREAK.join([
            'biom convert --to-tsv',
            f'-i {self.workdir}/feature-table.biom',
            f'-o {self.tsv}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)

    def remove_first_line(self):
        with open(self.tsv, 'r') as f:
            lines = f.readlines()
        with open(self.tsv, 'w') as f:
            f.writelines(lines[1:])


class ExportFeatureSequence(Export):

    feature_sequence_qza: str
    output_fa: str

    def main(self, feature_sequence_qza: str) -> str:

        self.feature_sequence_qza = feature_sequence_qza

        self.qza_to_fa()
        self.move_fa()

        return self.output_fa

    def qza_to_fa(self):
        self.qiime_tools_export(
            input_path=self.feature_sequence_qza,
            output_path=self.workdir)

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

    def main(self, taxonomy_qza: str) -> str:
        self.taxonomy_qza = taxonomy_qza

        self.qza_to_tsv()
        self.move_tsv()

        return self.tsv

    def qza_to_tsv(self):
        self.qiime_tools_export(
            input_path=self.taxonomy_qza,
            output_path=self.workdir)

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

    def main(self, aligned_sequence_qza: str) -> str:
        self.aligned_sequence_qza = aligned_sequence_qza

        self.qza_to_fa()
        self.move_fa()

        return self.output_fa

    def qza_to_fa(self):
        self.qiime_tools_export(
            input_path=self.aligned_sequence_qza,
            output_path=self.workdir)

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

    def main(self, tree_qza: str) -> str:
        self.tree_qza = tree_qza

        self.qza_to_nwk()
        self.move_nwk()

        return self.nwk

    def qza_to_nwk(self):
        self.qiime_tools_export(
            input_path=self.tree_qza,
            output_path=self.workdir)

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

    def main(self, distance_matrix_qza: str) -> str:
        self.distance_matrix_qza = distance_matrix_qza
        self.qza_to_tsv()
        self.move_tsv()
        return self.tsv

    def qza_to_tsv(self):
        self.qiime_tools_export(
            input_path=self.distance_matrix_qza,
            output_path=self.workdir)

    def move_tsv(self):
        self.tsv = edit_fpath(
            fpath=self.distance_matrix_qza,
            old_suffix='.qza',
            new_suffix='.tsv',
            dstdir=self.workdir
        )
        self.mv(f'{self.workdir}/distance-matrix.tsv', self.tsv)
