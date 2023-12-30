import os
import pandas as pd
from typing import Dict, Hashable
from .tools import edit_fpath
from .template import Processor
from .grouping import GROUP_COLUMN
from .lefse_plot_res import LefSePlotRes
from .lefse_plot_cladogram import LefSePlotCladogram


class LefSe(Processor):

    taxon_table_tsv_dict: Dict[str, str]
    sample_sheet: str
    colors: list

    def main(
            self,
            taxon_table_tsv_dict: Dict[str, str],
            sample_sheet: str,
            colors: list):

        self.taxon_table_tsv_dict = taxon_table_tsv_dict
        self.sample_sheet = sample_sheet
        self.colors = colors

        for level, tsv in self.taxon_table_tsv_dict.items():
            try:
                LefSeOneTaxonLevel(self.settings).main(
                    taxon_table_tsv=tsv,
                    taxon_level=level,
                    sample_sheet=self.sample_sheet,
                    colors=self.colors)
            except Exception as e:
                self.logger.warning(f'Failed to run LefSe on "{level}" taxon table, with Exception:\n{e}')


class LefSeOneTaxonLevel(Processor):

    DSTDIR_NAME = 'lefse'

    taxon_table_tsv: str
    taxon_level: str
    sample_sheet: str
    colors: list

    dstdir: str
    lefse_input: str
    lefse_result: str
    lefse_log: str
    feature_png: str
    cladogram_png: str

    def main(
            self,
            taxon_table_tsv: str,
            taxon_level: str,
            sample_sheet: str,
            colors: list):

        self.taxon_table_tsv = taxon_table_tsv
        self.taxon_level = taxon_level
        self.sample_sheet = sample_sheet
        self.colors = colors

        self.make_dstdir()
        self.add_taxon_level_prefix()
        self.insert_group_row()
        self.format_input()
        self.run_lefse()
        self.plot_feature()
        self.plot_cladogram()

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def add_taxon_level_prefix(self):
        self.taxon_table_tsv = AddTaxonLevelPrefix(self.settings).main(
            taxon_table_tsv=self.taxon_table_tsv)

    def insert_group_row(self):
        self.taxon_table_tsv = InsertGroupRow(self.settings).main(
            taxon_table_tsv=self.taxon_table_tsv,
            sample_sheet=self.sample_sheet)

    def format_input(self):
        self.lefse_input = f'{self.workdir}/lefse-{self.taxon_level}-input'
        self.lefse_log = f'{self.dstdir}/lefse-{self.taxon_level}.log'
        cmd = self.CMD_LINEBREAK.join([
            'lefse_format_input.py',
            self.taxon_table_tsv,
            self.lefse_input,
            f'-u 1',  # first row subject id
            f'-c 2',  # second row class id
            f'1>> {self.lefse_log}',
            f'2>> {self.lefse_log}'
        ])
        self.call(cmd)

    def run_lefse(self):
        self.lefse_result = f'{self.dstdir}/lefse-{self.taxon_level}-result.txt'
        cmd = self.CMD_LINEBREAK.join([
            'lefse_run.py',
            self.lefse_input,
            self.lefse_result,
            f'1>> {self.lefse_log}',
            f'2>> {self.lefse_log}'
        ])
        self.call(cmd)

    def plot_feature(self):
        self.feature_png = f'{self.dstdir}/lefse-{self.taxon_level}-features.png'
        LefSePlotRes().main(
            input_file=self.lefse_result,
            output_file=self.feature_png,
            colors=self.colors)

    def plot_cladogram(self):
        self.cladogram_png = f'{self.dstdir}/lefse-{self.taxon_level}-cladogram.png'

        LefSePlotCladogram().main(
            input_file=self.lefse_result,
            output_file=self.cladogram_png,
            colors=self.colors)


class AddTaxonLevelPrefix(Processor):

    LEVEL_PREFIXES = [
        'domain__',
        'phylum__',
        'class__',
        'order__',
        'family__',
        'genus__',
        'species__',
    ]

    taxon_table_tsv: str

    df: pd.DataFrame
    output_tsv: str

    def main(self, taxon_table_tsv: str) -> str:

        self.taxon_table_tsv = taxon_table_tsv

        self.read_taxon_table_tsv()
        self.df.index = map(self.add_level_prefixes, self.df.index)
        self.write_output_tsv()

        return self.output_tsv

    def add_level_prefixes(self, s: str) -> str:
        items = []
        for i, taxon in enumerate(s.split('|')):
            items.append(self.LEVEL_PREFIXES[i] + taxon)
        return '|'.join(items)

    def read_taxon_table_tsv(self):
        self.df = pd.read_csv(self.taxon_table_tsv, sep='\t', index_col=0)

    def write_output_tsv(self):
        self.output_tsv = edit_fpath(
            fpath=self.taxon_table_tsv,
            old_suffix='.tsv',
            new_suffix='-relabeled.tsv',
            dstdir=self.workdir
        )
        self.df.to_csv(self.output_tsv, sep='\t', index=True)


class InsertGroupRow(Processor):

    GROUP_INDEX = 'Group'
    NA_VALUE: str = 'None'

    taxon_table_tsv: str
    sample_sheet: str

    df: pd.DataFrame
    sample_to_group: Dict[Hashable, str]
    output_tsv: str

    def main(
            self,
            taxon_table_tsv: str,
            sample_sheet: str) -> str:

        self.taxon_table_tsv = taxon_table_tsv
        self.sample_sheet = sample_sheet

        self.read_taxon_table_tsv()
        self.set_sample_to_group()
        self.add_group_row()
        self.move_group_row()
        self.write_output_tsv()

        return self.output_tsv

    def read_taxon_table_tsv(self):
        self.df = pd.read_csv(self.taxon_table_tsv, sep='\t', index_col=0)

    def set_sample_to_group(self):
        sample_df = pd.read_csv(self.sample_sheet, index_col=0)
        self.sample_to_group = {}
        for sample, row in sample_df.iterrows():
            group = row[GROUP_COLUMN]
            if pd.isna(group):
                self.sample_to_group[sample] = self.NA_VALUE
            else:
                self.sample_to_group[sample] = group.replace(' ', '_')  # lefse does not allow space in group name

    def add_group_row(self):
        self.df.loc[self.GROUP_INDEX] = [
            self.sample_to_group[sample] for sample in self.df.columns
        ]

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
