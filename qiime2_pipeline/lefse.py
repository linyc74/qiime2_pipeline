import os
import pandas as pd
from typing import Dict, Hashable
from .utils import edit_fpath
from .template import Processor
from .grouping import GROUP_COLUMN
from .lefse_plot_res import LefSePlotRes
from .lefse_char_mapping import ORIGINAL_TO_NEW
from .lefse_plot_cladogram import LefSePlotCladogram


DSTDIR = 'lefse'


class LefSe(Processor):

    table_tsv_dict: Dict[str, str]
    sample_sheet: str
    colors: list

    def main(
            self,
            table_tsv_dict: Dict[str, str],
            sample_sheet: str,
            colors: list):

        self.table_tsv_dict = table_tsv_dict
        self.sample_sheet = sample_sheet
        self.colors = colors

        for name, tsv in self.table_tsv_dict.items():
            try:
                OneLefSe(self.settings).main(
                    table_tsv=tsv,
                    name=name,
                    sample_sheet=self.sample_sheet,
                    colors=self.colors)
            except Exception as e:
                self.logger.warning(f'Failed to run LefSe on "{name}" table, with Exception:\n{repr(e)}')


class OneLefSe(Processor):

    TAXON_LEVELS = [
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species'
    ]

    table_tsv: str
    name: str
    sample_sheet: str
    colors: list

    def main(
            self,
            table_tsv: str,
            name: str,
            sample_sheet: str,
            colors: list):

        self.table_tsv = table_tsv
        self.name = name
        self.sample_sheet = sample_sheet
        self.colors = colors

        os.makedirs(f'{self.outdir}/{DSTDIR}', exist_ok=True)

        self.lefse_features_workflow()
        if self.name in self.TAXON_LEVELS:  # do not run cladogram for the picrust2 pathway table
            self.lefse_cladogram_workflow()

    def lefse_features_workflow(self):
        table_tsv = FormatForLefse(self.settings).main(
            table_tsv=self.table_tsv,
            sample_sheet=self.sample_sheet,
            full_taxon_levels=False)

        lefse_input = f'{self.workdir}/lefse-{self.name}-features-input'
        self.lefse_format_input(
            table_tsv=table_tsv,
            lefse_input=lefse_input)

        lefse_result_tsv = f'{self.workdir}/lefse-{self.name}-features-result.tsv'
        self.lefse_run(
            lefse_input=lefse_input,
            lefse_result_tsv=lefse_result_tsv)

        # to avoid too many (unplottable) features
        limited_lefse_result_tsv = LimitLefseResult(self.settings).main(
            input_tsv=lefse_result_tsv)

        png = f'{self.outdir}/{DSTDIR}/lefse-{self.name}-features.png'
        LefSePlotRes().main(
            input_file=limited_lefse_result_tsv,
            output_file=png,
            sample_sheet=self.sample_sheet,
            colors=self.colors)

    def lefse_cladogram_workflow(self):
        table_tsv = FormatForLefse(self.settings).main(
            table_tsv=self.table_tsv,
            sample_sheet=self.sample_sheet,
            full_taxon_levels=True)

        lefse_input = f'{self.workdir}/lefse-{self.name}-cladogram-input'
        self.lefse_format_input(
            table_tsv=table_tsv,
            lefse_input=lefse_input)

        lefse_result_tsv = f'{self.workdir}/lefse-{self.name}-cladogram-result.tsv'
        self.lefse_run(
            lefse_input=lefse_input,
            lefse_result_tsv=lefse_result_tsv)

        png = f'{self.outdir}/{DSTDIR}/lefse-{self.name}-cladogram.png'
        LefSePlotCladogram().main(
            input_file=lefse_result_tsv,
            output_file=png,
            sample_sheet=self.sample_sheet,
            colors=self.colors)

    def lefse_format_input(self, table_tsv: str, lefse_input: str):
        log = f'{self.outdir}/lefse-{self.name}.log'
        cmd = self.CMD_LINEBREAK.join([
            'lefse_format_input.py',
            table_tsv,
            lefse_input,
            f'-u 1',  # first row subject id
            f'-c 2',  # second row class id
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)

    def lefse_run(self, lefse_input: str, lefse_result_tsv: str):
        log = f'{self.outdir}/lefse-{self.name}.log'
        cmd = self.CMD_LINEBREAK.join([
            'lefse_run.py',
            lefse_input,
            f'"{lefse_result_tsv}"',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)


class FormatForLefse(Processor):
    """
    When multiple taxon levels separated by '|', e.g. 'Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus|Streptococcus_mutans',
    LEfSe can compute LDA scores for all levels
    """

    LEVEL_PREFIXES = [
        'domain__',
        'phylum__',
        'class__',
        'order__',
        'family__',
        'genus__',
        'species__',
    ]
    GROUP_INDEX = 'Group'
    NA_VALUE: str = 'None'

    table_tsv: str
    sample_sheet: str
    full_taxon_levels: bool

    df: pd.DataFrame
    sample_to_group: Dict[Hashable, str]
    output_tsv: str

    def main(
            self,
            table_tsv: str,
            sample_sheet: str,
            full_taxon_levels: bool) -> str:

        self.table_tsv = table_tsv
        self.sample_sheet = sample_sheet
        self.full_taxon_levels = full_taxon_levels

        self.read_taxon_table_tsv()
        if self.full_taxon_levels:
            self.add_taxon_level_prefixes()
        else:
            self.shorten_taxon()
        self.set_sample_to_group()
        self.add_group_row()
        self.replace_invalid_characters()
        self.write_output_tsv()

        return self.output_tsv

    def read_taxon_table_tsv(self):
        self.df = pd.read_csv(self.table_tsv, sep='\t', index_col=0)

    def add_taxon_level_prefixes(self):

        def __add(s: str) -> str:
            items = []
            for i, taxon in enumerate(s.split('|')):
                items.append(self.LEVEL_PREFIXES[i] + taxon)
            return '|'.join(items)

        self.df.index = map(__add, self.df.index)

    def shorten_taxon(self):
        def shorten(s: str) -> str:
            return s.split('|')[-1]
        self.df.index = pd.Series(self.df.index).apply(shorten)

    def set_sample_to_group(self):
        sample_df = pd.read_csv(self.sample_sheet, index_col=0)
        self.sample_to_group = {}
        for sample, row in sample_df.iterrows():
            group = row[GROUP_COLUMN]
            if pd.isna(group):
                self.sample_to_group[sample] = self.NA_VALUE
            else:
                group = str(group)  # in case group is a number
                self.sample_to_group[sample] = group.replace(' ', '_')  # lefse does not allow space in group name

    def add_group_row(self):
        self.df.loc[self.GROUP_INDEX] = [
            self.sample_to_group[sample] for sample in self.df.columns
        ]

        # move to the second row
        indexes = list(self.df.index)
        indexes = [indexes[-1]] + indexes[:-1]
        self.df = self.df.loc[indexes]

    def replace_invalid_characters(self):
        def replace(s: str) -> str:
            for a, b in ORIGINAL_TO_NEW:
                s = s.replace(a, b)
            return s
        self.df.index = pd.Series(self.df.index).apply(replace)

    def write_output_tsv(self):
        suffix = 'full-taxon' if self.full_taxon_levels else 'short-taxon'
        self.output_tsv = edit_fpath(
            fpath=self.table_tsv,
            old_suffix='.tsv',
            new_suffix=f'-lefse-{suffix}.tsv',
            dstdir=self.workdir
        )
        self.df.to_csv(self.output_tsv, sep='\t', index=True)


class LimitLefseResult(Processor):

    MAX_FEATURES_PER_GROUP = 50

    input_tsv: str

    indf: pd.DataFrame
    outdf: pd.DataFrame
    output_tsv: str

    def main(self, input_tsv: str) -> str:
        self.input_tsv = input_tsv

        self.init_indf_outdf()

        for group in self.indf['Group'].unique():
            if pd.isna(group):
                continue
            self.keep_the_highest_lda_for(group)

        self.save_output_tsv()

        return self.output_tsv

    def init_indf_outdf(self):
        self.indf = pd.read_csv(
            self.input_tsv,
            sep='\t',
            header=None,
            names=['Column 1', 'Group', 'LDA', 'Column 4'],
            index_col=0
        )
        self.outdf = pd.DataFrame(columns=self.indf.columns)

    def keep_the_highest_lda_for(self, group: str):
        subdf = self.indf[self.indf['Group'] == group]
        subdf = subdf.sort_values(
            by='LDA',
            ascending=False
        ).head(
            n=self.MAX_FEATURES_PER_GROUP
        )
        self.outdf = pd.concat([self.outdf, subdf])

    def save_output_tsv(self):
        self.output_tsv = edit_fpath(
            fpath=self.input_tsv,
            old_suffix='.tsv',
            new_suffix='-filtered.tsv',
            dstdir=self.workdir)

        self.outdf.to_csv(self.output_tsv, sep='\t', index=True, header=False)
