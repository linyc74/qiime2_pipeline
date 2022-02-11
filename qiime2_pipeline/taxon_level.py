import os
import pandas as pd
from typing import Dict
from .template import Processor, Settings


class TaxonLevel(Processor):

    OUTPUT_DIRNAME = 'taxon-level'
    TAXON_KEY_TO_LEVEL = {
        'phylum': 2,
        'class': 3,
        'order': 4,
        'family': 5,
        'genus': 6,
        'species': 7
    }

    df: pd.DataFrame

    collapsed_dfs: Dict[str, pd.DataFrame]
    tsv_dict: Dict[str, str]

    def __init__(self, settings: Settings):
        super().__init__(settings)

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
        for key, level in self.TAXON_KEY_TO_LEVEL.items():
            self.collapsed_dfs[key] = self.get_collapsed_df(level=level)

    def get_collapsed_df(self, level: int) -> pd.DataFrame:
        df = self.df.copy(deep=True)
        df.index = [
            trim_taxon(taxon=idx, level=level)
            for idx in self.df.index
        ]
        df = df.groupby(level=0).sum()  # collapse by 0-level index
        return df

    def save_output_tsvs(self):
        dstdir = f'{self.outdir}/{self.OUTPUT_DIRNAME}'
        os.makedirs(dstdir, exist_ok=True)

        self.tsv_dict = {}
        for key, df in self.collapsed_dfs.items():
            tsv = f'{dstdir}/{key}.tsv'
            df.to_csv(tsv, sep='\t', index=True)
            self.tsv_dict[key] = tsv


def feature_label_to_taxon(s: str) -> str:
    """
    input: 'X; x__AAA; x__BBB; x__CCC; x__DDD; x__EEE; x__FFF; x__GGG'
    output: 'AAA|BBB|CCC|DDD|EEE|FFF|GGG'
    """
    items = s.split('; ')[1:]
    return '|'.join(i[3:] for i in items)


def trim_taxon(taxon: str, level: int) -> str:
    """
    taxon_hierarchy: '1|2|3|4|5|6|7'
    level: 5
    output: '1|2|3|4|5'
    """
    return '|'.join(taxon.split('|')[:level])

