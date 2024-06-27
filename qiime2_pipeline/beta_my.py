import pandas as pd
from abc import ABC
from .utils import edit_fpath
from .template import Processor
from .normalization import CountNormalization
from .embedding_core_process import PCACore, TSNECore, EmbeddingProcessTemplate


class MyBetaDiversity(Processor):

    feature_table_tsv: str
    sample_sheet: str
    colors: list

    def main(
            self,
            feature_table_tsv: str,
            sample_sheet: str,
            colors: list):

        self.feature_table_tsv = feature_table_tsv
        self.sample_sheet = sample_sheet
        self.colors = colors

        for Process in [PCAProcess, TSNEProcess]:
            Process(self.settings).main(
                tsv=self.feature_table_tsv,
                sample_sheet=self.sample_sheet,
                colors=self.colors)


class EmbeddingProcess(EmbeddingProcessTemplate, ABC):

    DSTDIR_NAME = 'feature-embedding'
    LOG_PSEUDOCOUNT = True
    NORMALIZE_BY_SAMPLE_READS = False

    def preprocessing(self):
        self.df = CountNormalization(self.settings).main(
            df=self.df,
            log_pseudocount=self.LOG_PSEUDOCOUNT,
            by_sample_reads=self.NORMALIZE_BY_SAMPLE_READS
        )


class PCAProcess(EmbeddingProcess):

    NAME = 'PCA'
    XY_COLUMNS = ['PC 1', 'PC 2']

    proportion_explained_series: pd.Series

    def main(
            self,
            tsv: str,
            sample_sheet: str,
            colors: list):

        self.tsv = tsv
        self.sample_sheet = sample_sheet
        self.colors = colors

        self.run_main_workflow()
        self.write_proportion_explained()

    def embedding(self):
        self.sample_coordinate_df, self.proportion_explained_series = PCACore(self.settings).main(
            df=self.df,
            data_structure='row_features'
        )

    def write_proportion_explained(self):
        tsv = edit_fpath(
            fpath=self.tsv,
            old_suffix='.tsv',
            new_suffix=f'-{self.NAME.lower()}-proportion-explained.tsv',
            dstdir=self.dstdir
        )
        self.proportion_explained_series.to_csv(
            tsv,
            sep='\t',
            header=['Proportion Explained']
        )


class TSNEProcess(EmbeddingProcess):

    NAME = 't-SNE'
    XY_COLUMNS = ['t-SNE 1', 't-SNE 2']

    def main(
            self,
            tsv: str,
            sample_sheet: str,
            colors: list):

        self.tsv = tsv
        self.sample_sheet = sample_sheet
        self.colors = colors

        self.run_main_workflow()

    def embedding(self):
        self.sample_coordinate_df = TSNECore(self.settings).main(
            df=self.df,
            data_structure='row_features'
        )

