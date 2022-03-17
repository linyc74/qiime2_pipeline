import numpy as np
import pandas as pd
from .template import Processor


class CountNormalization(Processor):

    df: pd.DataFrame
    by_sample_reads: bool
    sample_reads_unit: int
    log_pseudocount: bool

    def main(
            self,
            df: pd.DataFrame,
            by_sample_reads: bool,
            sample_reads_unit: int,
            log_pseudocount: bool) -> pd.DataFrame:

        self.df = df
        self.by_sample_reads = by_sample_reads
        self.sample_reads_unit = sample_reads_unit
        self.log_pseudocount = log_pseudocount

        self.pseudocount()
        self.normalize_by_sample_reads()
        self.log10()

        return self.df

    def pseudocount(self):
        if self.log_pseudocount:
            self.df = self.df + 1

    def normalize_by_sample_reads(self):
        if self.by_sample_reads:
            sum_per_column = np.sum(self.df, axis=0) / self.sample_reads_unit
            self.df = np.divide(self.df, sum_per_column)

    def log10(self):
        if self.log_pseudocount:
            self.df = np.log10(self.df)
