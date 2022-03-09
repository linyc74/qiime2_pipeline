import os
import pandas as pd
from typing import List, Dict
from .tools import edit_fpath
from .template import Processor, Settings


class LefSe(Processor):

    taxon_table_tsv_dict: Dict[str, str]
    group_keywords: List[str]

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            taxon_table_tsv_dict: Dict[str, str],
            group_keywords: List[str]):

        self.taxon_table_tsv_dict = taxon_table_tsv_dict
        self.group_keywords = group_keywords

        for level, tsv in self.taxon_table_tsv_dict.items():
            LefSeOneTaxonLevel(self.settings).main(
                taxon_table_tsv=tsv,
                taxon_level=level,
                group_keywords=self.group_keywords)


class LefSeOneTaxonLevel(Processor):

    DSTDIR_NAME = 'lefse'
    DPI = 600
    FIGURE_WIDTH = 8

    taxon_table_tsv: str
    taxon_level: str
    group_keywords: List[str]

    dstdir: str
    lefse_input: str
    lefse_result: str
    lefse_log: str
    feature_png: str
    cladogram_png: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            taxon_table_tsv: str,
            taxon_level: str,
            group_keywords: List[str]):

        self.taxon_table_tsv = taxon_table_tsv
        self.taxon_level = taxon_level
        self.group_keywords = group_keywords

        self.make_dstdir()
        self.insert_group_row()
        self.format_input()
        self.run_lefse()
        self.plot_feature()
        self.plot_cladogram()

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def insert_group_row(self):
        self.taxon_table_tsv = InsertGroupRow(self.settings).main(
            taxon_table_tsv=self.taxon_table_tsv,
            group_keywords=self.group_keywords)

    def format_input(self):
        self.lefse_input = f'{self.workdir}/lefse-{self.taxon_level}-input'
        self.lefse_log = f'{self.dstdir}/lefse-{self.taxon_level}.log'
        cmd = self.CMD_LINEBREAK.join([
            'lefse_format_input.py',
            self.taxon_table_tsv,
            self.lefse_input,
            f'-u 1',  # first row subject id
            f'-c 2',  # second row class id
            f'1> {self.lefse_log} 2> {self.lefse_log}'
        ])
        self.call(cmd)

    def run_lefse(self):
        self.lefse_result = f'{self.dstdir}/lefse-{self.taxon_level}-result.txt'
        cmd = self.CMD_LINEBREAK.join([
            'lefse_run.py',
            self.lefse_input,
            self.lefse_result,
            f'1> {self.lefse_log} 2> {self.lefse_log}'
        ])
        self.call(cmd)

    def plot_feature(self):
        self.feature_png = f'{self.dstdir}/lefse-{self.taxon_level}-features.png'
        cmd = self.CMD_LINEBREAK.join([
            'lefse_plot_res.py',
            f'--dpi {self.DPI}',
            f'--width {self.FIGURE_WIDTH}',
            self.lefse_result,
            self.feature_png,
            f'1> {self.lefse_log} 2> {self.lefse_log}'
        ])
        self.call(cmd)

    def plot_cladogram(self):
        self.cladogram_png = f'{self.dstdir}/lefse-{self.taxon_level}-cladogram.png'
        cmd = self.CMD_LINEBREAK.join([
            'lefse_plot_cladogram.py',
            f'--dpi {self.DPI}',
            '--format png',
            self.lefse_result,
            self.cladogram_png,
            f'1> {self.lefse_log} 2> {self.lefse_log}'
        ])
        self.call(cmd)


class InsertGroupRow(Processor):

    GROUP_INDEX = 'Group'
    NA_VALUE: str = 'None'

    taxon_table_tsv: str
    group_keywords: List[str]

    df: pd.DataFrame
    output_tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            taxon_table_tsv: str,
            group_keywords: List[str]) -> str:

        self.taxon_table_tsv = taxon_table_tsv
        self.group_keywords = group_keywords

        self.read_taxon_table_tsv()
        self.add_group_row()
        self.move_group_row()
        self.write_output_tsv()

        return self.output_tsv

    def read_taxon_table_tsv(self):
        self.df = pd.read_csv(self.taxon_table_tsv, sep='\t', index_col=0)

    def add_group_row(self):
        self.df.loc[self.GROUP_INDEX] = [
            self.__to_group(c) for c in self.df.columns
        ]

    def __to_group(self, sample_name: str) -> str:
        for k in self.group_keywords:
            if k in sample_name:
                return k
        return self.NA_VALUE

    def move_group_row(self):
        indexes = list(self.df.index)
        indexes = [indexes[-1]] + indexes[:-1]
        self.df = self.df.loc[indexes]

    def write_output_tsv(self):
        self.output_tsv = edit_fpath(
            fpath=self.taxon_table_tsv,
            old_suffix='.tsv',
            new_suffix='-grouped.tsv',
            dstdir=self.workdir
        )
        self.df.to_csv(self.output_tsv, sep='\t', index=True)
