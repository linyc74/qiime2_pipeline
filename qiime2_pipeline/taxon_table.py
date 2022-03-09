import os
import pandas as pd
from typing import Dict
from .template import Processor, Settings


class TaxonTable(Processor):

    DSTDIR_NAME = 'taxon-table'
    TAXON_LEVELS = [
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
    ]

    df: pd.DataFrame

    collapsed_dfs: Dict[str, pd.DataFrame]
    tsv_dict: Dict[str, str]

    def main(self, labeled_feature_table_tsv: str) -> Dict[str, str]:
        self.df = pd.read_csv(labeled_feature_table_tsv, sep='\t', index_col=0)

        self.feature_label_to_taxon()
        self.collapse_taxon_at_various_levels()
        self.save_output_tsvs()

        return self.tsv_dict

    def feature_label_to_taxon(self):
        self.df.index = [feature_label_to_taxon(idx) for idx in self.df.index]

    def collapse_taxon_at_various_levels(self):
        self.collapsed_dfs = {}
        for level in self.TAXON_LEVELS:
            self.collapsed_dfs[level] = CollapseTaxon(self.settings).main(
                df=self.df,
                level=level)

    def save_output_tsvs(self):
        dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(dstdir, exist_ok=True)

        self.tsv_dict = {}
        for level, df in self.collapsed_dfs.items():
            tsv = f'{dstdir}/{level}-table.tsv'
            df.to_csv(tsv, sep='\t', index=True)
            self.tsv_dict[level] = tsv


class CollapseTaxon(Processor):

    TAXON_LEVEL_STR_TO_INT = {
        'phylum': 2,
        'class': 3,
        'order': 4,
        'family': 5,
        'genus': 6,
        'species': 7
    }

    df: pd.DataFrame
    level: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, df: pd.DataFrame, level: str) -> pd.DataFrame:
        self.df = df.copy(deep=True)
        self.level = level

        self.trim_taxon()
        self.collapse_by_index()

        return self.df

    def trim_taxon(self):
        int_level = self.TAXON_LEVEL_STR_TO_INT[self.level]
        self.df.index = [
            trim_taxon(taxon=taxon, int_level=int_level)
            for taxon in self.df.index
        ]

    def collapse_by_index(self):
        # groupby automatically sorts the index
        self.df = self.df.groupby(level=0).sum()  # 0-level index


def feature_label_to_taxon(s: str) -> str:
    """
    input: 'X; x__AAA; x__BBB; x__CCC; x__DDD; x__EEE; x__FFF; x__GGG'
    output: 'AAA|BBB|CCC|DDD|EEE|FFF|GGG'
    """
    items = s.split('; ')[1:]
    return '|'.join(i[3:] for i in items)


def trim_taxon(taxon: str, int_level: int) -> str:
    """
    taxon_hierarchy: '1|2|3|4|5|6|7'
    level: 5
    output: '1|2|3|4|5'
    """
    return '|'.join(taxon.split('|')[:int_level])
