import os
import pandas as pd
from typing import Dict
from .utils import edit_fpath
from .template import Processor
from .fasta import FastaParser, FastaWriter


class PICRUSt2(Processor):

    DSTDIR_NAME = 'picrust2'

    labeled_feature_sequence_fa: str
    labeled_feature_table_tsv: str

    feature_sequence_fa: str
    feature_table_tsv: str
    picrust2_workdir: str
    tsv_dict: Dict[str, str]

    def main(self, labeled_feature_sequence_fa: str, labeled_feature_table_tsv: str) -> Dict[str, str]:

        self.labeled_feature_sequence_fa = labeled_feature_sequence_fa
        self.labeled_feature_table_tsv = labeled_feature_table_tsv

        try:
            self.remove_white_space_in_headers()
            self.remove_existing_workdir()
            self.run_picrust2()
            self.unzip_and_add_descriptions()
            return self.tsv_dict
        except Exception as e:
            self.logger.warning(e)
            self.clean_up_if_failed()
            return {}

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

    def remove_existing_workdir(self):
        self.picrust2_workdir = f'{self.workdir}/{self.DSTDIR_NAME}'
        if os.path.exists(self.picrust2_workdir):
            self.call(f'rm -r {self.picrust2_workdir}')  # to avoid directory already exists error

    def run_picrust2(self):
        log = f'{self.outdir}/picrust2.log'
        args = [
            'picrust2_pipeline.py',
            f'--study_fasta {self.feature_sequence_fa}',
            f'--input {self.feature_table_tsv}',
            f'--output {self.picrust2_workdir}',
            f'--processes {self.threads}',
            f'1>> "{log}"',
            f'2>> "{log}"',
        ]
        self.call(self.CMD_LINEBREAK.join(args))

    def unzip_and_add_descriptions(self):
        tsv_dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(tsv_dstdir, exist_ok=True)
        self.tsv_dict = {}

        for key, in_tsv, map_type in [
            (
                'picrust2-pathway',
                f'{self.picrust2_workdir}/pathways_out/path_abun_unstrat.tsv.gz',
                'METACYC',
            ),
            (
                'picrust2-EC',
                f'{self.picrust2_workdir}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz',
                'EC',
            ),
            (
                'picrust2-KEGG-ortholog',
                f'{self.picrust2_workdir}/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz',
                'KO',
            ),
        ]:
            out_tsv = f'{tsv_dstdir}/{key}-table.tsv'
            self.add_descriptions(in_tsv, map_type, out_tsv)
            self.tsv_dict[key] = out_tsv

    def add_descriptions(self, in_tsv: str, map_type: str, out_tsv: str):
        log = f'{self.outdir}/picrust2.log'
        args = [
            'add_descriptions.py',
            f'--input "{in_tsv}"',
            f'--map_type {map_type}',
            f'--output "{out_tsv}"',
            f'1>> "{log}"',
            f'2>> "{log}"',
        ]
        self.call(self.CMD_LINEBREAK.join(args))

        df = pd.read_csv(out_tsv, sep='\t', index_col=0)
        df = merge_description_into_index(df=df)
        df.to_csv(out_tsv, sep='\t')

    def clean_up_if_failed(self):
        d = f'{self.outdir}/{self.DSTDIR_NAME}'
        if os.path.exists(d):
            self.call(f'rm -r {self.outdir}/{self.DSTDIR_NAME}')


def merge_description_into_index(df: pd.DataFrame) -> pd.DataFrame:
    index_name = df.index.name
    df.index = [f'[{idx}]' for idx in df.index]
    df.index.name = index_name
    df.index = df.index.str.cat(df['description'], sep=' ')
    return df.drop(columns='description')
