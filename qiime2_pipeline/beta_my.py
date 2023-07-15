import pandas as pd
from abc import ABC
from typing import List
from .tools import edit_fpath
from .template import Processor
from .normalization import CountNormalization
from .embedding_core_process import PCACore, NMDSCore, TSNECore, EmbeddingProcessTemplate


class MyBetaDiversity(Processor):

    feature_table_tsv: str
    group_keywords: List[str]

    def main(
            self,
            feature_table_tsv: str,
            group_keywords: List[str]):

        self.feature_table_tsv = feature_table_tsv
        self.group_keywords = group_keywords

        for Process in [PCAProcess, NMDSProcess, TSNEProcess]:
            Process(self.settings).main(
                tsv=self.feature_table_tsv,
                group_keywords=self.group_keywords)


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
    XY_COLUMNS = ('PC 1', 'PC 2')

    proportion_explained_series: pd.Series

    def main(
            self,
            tsv: str,
            group_keywords: List[str]):

        self.tsv = tsv
        self.group_keywords = group_keywords

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


class NMDSProcess(EmbeddingProcess):

    NAME = 'NMDS'
    XY_COLUMNS = ('NMDS 1', 'NMDS 2')

    stress: float

    def main(
            self,
            tsv: str,
            group_keywords: List[str]):

        self.tsv = tsv
        self.group_keywords = group_keywords

        self.run_main_workflow()
        self.write_stress()

    def embedding(self):
        self.sample_coordinate_df, self.stress = NMDSCore(self.settings).main(
            df=self.df,
            data_structure='row_features'
        )

    def write_stress(self):
        txt = edit_fpath(
            fpath=self.tsv,
            old_suffix='.tsv',
            new_suffix='-nmds-stress.txt',
            dstdir=self.dstdir)
        with open(txt, 'w') as fh:
            fh.write(str(self.stress))


class TSNEProcess(EmbeddingProcess):

    NAME = 't-SNE'
    XY_COLUMNS = ('t-SNE 1', 't-SNE 2')

    def main(
            self,
            tsv: str,
            group_keywords: List[str]):

        self.tsv = tsv
        self.group_keywords = group_keywords

        self.run_main_workflow()

    def embedding(self):
        self.sample_coordinate_df = TSNECore(self.settings).main(
            df=self.df,
            data_structure='row_features'
        )

