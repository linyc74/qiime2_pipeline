import pandas as pd
from os.path import exists
from .tools import edit_fpath
from .template import Processor
from .fasta import FastaParser, FastaWriter


class PICRUSt2(Processor):

    DSTDIR_NAME = 'picrust2'

    labeled_feature_sequence_fa: str
    labeled_feature_table_tsv: str

    feature_sequence_fa: str
    feature_table_tsv: str
    pathway_table_tsv: str

    def main(self, labeled_feature_sequence_fa: str, labeled_feature_table_tsv: str) -> str:

        self.labeled_feature_sequence_fa = labeled_feature_sequence_fa
        self.labeled_feature_table_tsv = labeled_feature_table_tsv

        self.remove_white_space_in_headers()
        self.remove_existing_dstdir()
        self.run_picrust2()
        self.unzip_pathway_table_tsv()
        self.move_picrust2_dir()

        return self.pathway_table_tsv

    def remove_white_space_in_headers(self):

        self.feature_sequence_fa = edit_fpath(
            fpath=self.labeled_feature_sequence_fa,
            old_suffix='.fa',
            new_suffix='-remove-white-space.fa',
            dstdir=self.workdir
        )

        with FastaParser(self.labeled_feature_sequence_fa) as parser:
            with FastaWriter(self.feature_sequence_fa) as writer:
                for header, sequence in parser:
                    writer.write(header.split(';')[0], sequence)

        self.feature_table_tsv = edit_fpath(
            fpath=self.labeled_feature_table_tsv,
            old_suffix='.tsv',
            new_suffix='-remove-white-space.tsv',
            dstdir=self.workdir
        )

        df = pd.read_csv(self.labeled_feature_table_tsv, sep='\t', index_col=0)
        df.index = df.index.str.split(';').str[0]
        df.to_csv(self.feature_table_tsv, sep='\t')

    def remove_existing_dstdir(self):
        if exists(f'{self.workdir}/{self.DSTDIR_NAME}'):
            self.call(f'rm -r {self.workdir}/{self.DSTDIR_NAME}')  # to avoid directory already exists error

    def run_picrust2(self):
        log = f'{self.outdir}/picrust2.log'
        args = [
            'picrust2_pipeline.py',
            f'--study_fasta {self.feature_sequence_fa}',
            f'--input {self.feature_table_tsv}',
            f'--output {self.workdir}/{self.DSTDIR_NAME}',
            f'--processes {self.threads}',
            f'1>> "{log}"',
            f'2>> "{log}"',
        ]
        self.call(self.CMD_LINEBREAK.join(args))

    def unzip_pathway_table_tsv(self):
        self.pathway_table_tsv = f'{self.outdir}/picrust2-pathway-table.tsv'
        args = [
            'gunzip --keep --stdout',
            f'{self.workdir}/{self.DSTDIR_NAME}/pathways_out/path_abun_unstrat.tsv.gz',
            f'> "{self.pathway_table_tsv}"'
        ]
        self.call(self.CMD_LINEBREAK.join(args))

    def move_picrust2_dir(self):
        self.call(f'mv "{self.workdir}/{self.DSTDIR_NAME}" "{self.outdir}/"')
