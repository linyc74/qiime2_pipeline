import pandas as pd
from typing import List
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from .tools import edit_fpath
from .template import Processor, Settings
from .embedding_core import NMDSCore, TSNECore
from .embedding_process import EmbeddingProcessTemplate


class PCoAProcess(EmbeddingProcessTemplate):

    NAME = 'PCoA'
    XY_COLUMNS = ('PC1', 'PC2')
    DSTDIR_NAME = 'beta-embedding'

    proportion_explained_serise: pd.Series

    def main(
            self,
            tsv: str,
            group_keywords: List[str]):

        self.tsv = tsv
        self.group_keywords = group_keywords

        self.run_main_workflow()
        self.write_proportion_explained()

    def run_embedding(self):
        df = self.distance_matrix
        dist_mat = DistanceMatrix(df, list(df.columns))
        result = pcoa(distance_matrix=dist_mat)
        self.sample_coordinate_df = result.samples
        self.proportion_explained_serise = result.proportion_explained

    def write_proportion_explained(self):
        tsv = edit_fpath(
            fpath=self.tsv,
            old_suffix='.tsv',
            new_suffix=f'-{self.NAME.lower()}-proportion-explained.tsv',
            dstdir=self.dstdir
        )
        self.proportion_explained_serise.to_csv(
            tsv,
            sep='\t',
            header=['Proportion Explained']
        )


class NMDSProcess(EmbeddingProcessTemplate):

    NAME = 'NMDS'
    XY_COLUMNS = ('NMDS 1', 'NMDS 2')
    DSTDIR_NAME = 'beta-embedding'

    stress: float

    def main(
            self,
            tsv: str,
            group_keywords: List[str]):

        self.tsv = tsv
        self.group_keywords = group_keywords

        self.run_main_workflow()
        self.write_stress()

    def run_embedding(self):
        self.sample_coordinate_df, self.stress = NMDSCore(self.settings).main(
            df=self.distance_matrix,
            data_structure='distance_matrix'
        )

    def write_stress(self):
        txt = edit_fpath(
            fpath=self.tsv,
            old_suffix='.tsv',
            new_suffix='-nmds-stress.txt',
            dstdir=self.dstdir)
        with open(txt, 'w') as fh:
            fh.write(str(self.stress))


class TSNEProcess(EmbeddingProcessTemplate):

    NAME = 't-SNE'
    XY_COLUMNS = ('t-SNE 1', 't-SNE 2')
    DSTDIR_NAME = 'beta-embedding'

    def main(
            self,
            tsv: str,
            group_keywords: List[str]):

        self.tsv = tsv
        self.group_keywords = group_keywords
        self.run_main_workflow()

    def run_embedding(self):
        self.sample_coordinate_df = TSNECore(self.settings).main(
            df=self.distance_matrix,
            data_structure='distance_matrix'
        )


class BatchEmbeddingProcess(Processor):

    distance_matrix_tsvs: List[str]
    group_keywords: List[str]

    embedding: EmbeddingProcessTemplate

    def main(
            self,
            distance_matrix_tsvs: List[str],
            group_keywords: List[str]):

        self.distance_matrix_tsvs = distance_matrix_tsvs
        self.group_keywords = group_keywords

        for tsv in self.distance_matrix_tsvs:
            self.logger.debug(f'{self.embedding.NAME} for {tsv}')
            self.embedding.main(
                tsv=tsv,
                group_keywords=self.group_keywords)


class BatchPCoAProcess(BatchEmbeddingProcess):

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.embedding = PCoAProcess(self.settings)


class BatchNMDSProcess(BatchEmbeddingProcess):

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.embedding = NMDSProcess(self.settings)


class BatchTSNEProcess(BatchEmbeddingProcess):

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.embedding = TSNEProcess(self.settings)
