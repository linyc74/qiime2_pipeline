import numpy as np
import pandas as pd
from typing import Tuple
from .tools import edit_fpath
from .template import Processor
from .importing import ImportFeatureTable
from .exporting import ExportFeatureTable


class FeatureNormalization(Processor):

    feature_table_qza: str
    log_pseudocount: bool

    feature_table_tsv: str
    df: pd.DataFrame
    normalized_feature_table_qza: str
    normalized_feature_table_tsv: str

    def main(
            self,
            feature_table_qza: str,
            log_pseudocount: bool) -> Tuple[str, str]:

        self.feature_table_qza = feature_table_qza
        self.log_pseudocount = log_pseudocount

        self.decompress()
        self.read_tsv()
        if log_pseudocount:
            self.df = log10_pseudocount(self.df)
        self.save_tsv()
        self.compress()

        return self.normalized_feature_table_tsv, self.normalized_feature_table_qza

    def decompress(self):
        self.feature_table_tsv = ExportFeatureTable(self.settings).main(
            feature_table_qza=self.feature_table_qza)

    def read_tsv(self):
        self.df = pd.read_csv(
            self.feature_table_tsv,
            sep='\t',
            index_col=0,
            skiprows=1  # exclude 1st line from qza (# Constructed from biom file)
        )

    def save_tsv(self):
        self.normalized_feature_table_tsv = edit_fpath(
            fpath=self.feature_table_qza,
            old_suffix='.qza',
            new_suffix='-normalized.tsv',
            dstdir=self.outdir)
        self.df.to_csv(self.normalized_feature_table_tsv, sep='\t', index=True)

    def compress(self):
        self.normalized_feature_table_qza = ImportFeatureTable(self.settings).main(
            feature_table_tsv=self.normalized_feature_table_tsv)


def log10_pseudocount(df: pd.DataFrame) -> pd.DataFrame:
    return np.log10(df + 1)
