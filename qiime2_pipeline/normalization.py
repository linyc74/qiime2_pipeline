import numpy as np
import pandas as pd
from .template import Processor


class CountNormalization(Processor):

    df: pd.DataFrame
    log_pseudocount: bool
    by_sample_reads: bool
    sample_reads_unit: int

    def main(
            self,
            df: pd.DataFrame,
            log_pseudocount: bool,
            by_sample_reads: bool,
            sample_reads_unit: int = 10000) -> pd.DataFrame:

        self.df = df
        self.log_pseudocount = log_pseudocount
        self.by_sample_reads = by_sample_reads
        self.sample_reads_unit = sample_reads_unit

        self.normalize_by_sample_reads()
        self.pseudocount_then_log10()

        return self.df

    def normalize_by_sample_reads(self):
        if self.by_sample_reads:
            sum_per_column = np.sum(self.df, axis=0) / self.sample_reads_unit
            self.df = np.divide(self.df, sum_per_column)

    def pseudocount_then_log10(self):
        if self.log_pseudocount:
            self.df = np.log10(self.df + 1)
