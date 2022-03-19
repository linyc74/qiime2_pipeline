import os
import pandas as pd
from typing import List, Tuple
from skbio import DistanceMatrix
from abc import ABC, abstractmethod
from skbio.stats.ordination import pcoa
from .tools import edit_fpath
from .grouping import AddGroupColumn
from .embedding_core import ScatterPlot
from .template import Processor, Settings
from .embedding_core import NMDSCore, TSNECore


class BetaEmbedding(Processor, ABC):

    NAME: str
    XY_COLUMNS: Tuple[str, str]
    GROUP_COLUMN: str = AddGroupColumn.GROUP_COLUMN

    distance_matrix_tsv: str
    group_keywords: List[str]

    distance_matrix: pd.DataFrame
    sample_coordinate_df: pd.DataFrame
    dstdir: str

    @abstractmethod
    def main(
            self,
            distance_matrix_tsv: str,
            group_keywords: List[str]):
        pass

    def run_main_workflow(self):
        self.load_distance_matrix()
        self.run_embedding()
        self.add_group_column()
        self.make_dstdir()
        self.write_sample_coordinate()
        self.plot_sample_coordinate()

    def load_distance_matrix(self):
        self.distance_matrix = pd.read_csv(
            self.distance_matrix_tsv,
            sep='\t',
            index_col=0)

    @abstractmethod
    def run_embedding(self):
        pass

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/beta-embedding'
        os.makedirs(self.dstdir, exist_ok=True)

    def add_group_column(self):
        self.sample_coordinate_df = AddGroupColumn(self.settings).main(
            df=self.sample_coordinate_df,
            group_keywords=self.group_keywords)

    def write_sample_coordinate(self):
        tsv = self.__get_sample_coordinate_fpath(ext='tsv')
        self.sample_coordinate_df.to_csv(tsv, sep='\t')

    def plot_sample_coordinate(self):
        png = self.__get_sample_coordinate_fpath(ext='png')
        ScatterPlot(self.settings).main(
            sample_coordinate_df=self.sample_coordinate_df,
            x_column=self.XY_COLUMNS[0],
            y_column=self.XY_COLUMNS[1],
            hue_column=self.GROUP_COLUMN,
            output_png=png)

    def __get_sample_coordinate_fpath(self, ext: str) -> str:
        name = self.NAME.lower().replace('-', '')
        return edit_fpath(
            fpath=self.distance_matrix_tsv,
            old_suffix='.tsv',
            new_suffix=f'-{name}-sample-coordinate.{ext}',
            dstdir=self.dstdir
        )


class PCoA(BetaEmbedding):

    NAME = 'PCoA'
    XY_COLUMNS = ('PC1', 'PC2')

    proportion_explained_serise: pd.Series

    def main(
            self,
            distance_matrix_tsv: str,
            group_keywords: List[str]):

        self.distance_matrix_tsv = distance_matrix_tsv
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
            fpath=self.distance_matrix_tsv,
            old_suffix='.tsv',
            new_suffix=f'-{self.NAME.lower()}-proportion-explained.tsv',
            dstdir=self.dstdir
        )
        self.proportion_explained_serise.to_csv(
            tsv,
            sep='\t',
            header=['Proportion Explained']
        )


class NMDS(BetaEmbedding):

    NAME = 'NMDS'
    XY_COLUMNS = ('NMDS 1', 'NMDS 2')

    stress: float

    def main(
            self,
            distance_matrix_tsv: str,
            group_keywords: List[str]):

        self.distance_matrix_tsv = distance_matrix_tsv
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
            fpath=self.distance_matrix_tsv,
            old_suffix='.tsv',
            new_suffix='-nmds-stress.txt',
            dstdir=self.dstdir)
        with open(txt, 'w') as fh:
            fh.write(str(self.stress))


class TSNE(BetaEmbedding):

    NAME = 't-SNE'
    XY_COLUMNS = ('t-SNE 1', 't-SNE 2')

    def main(
            self,
            distance_matrix_tsv: str,
            group_keywords: List[str]):

        self.distance_matrix_tsv = distance_matrix_tsv
        self.group_keywords = group_keywords
        self.run_main_workflow()

    def run_embedding(self):
        self.sample_coordinate_df = TSNECore(self.settings).main(
            df=self.distance_matrix,
            data_structure='distance_matrix'
        )


class BatchBetaEmbedding(Processor):

    distance_matrix_tsvs: List[str]
    group_keywords: List[str]

    beta_embedding: BetaEmbedding

    def main(
            self,
            distance_matrix_tsvs: List[str],
            group_keywords: List[str]):

        self.distance_matrix_tsvs = distance_matrix_tsvs
        self.group_keywords = group_keywords

        for tsv in self.distance_matrix_tsvs:
            self.logger.debug(f'{self.beta_embedding.NAME} for {tsv}')
            self.beta_embedding.main(
                distance_matrix_tsv=tsv,
                group_keywords=self.group_keywords)


class BatchPCoA(BatchBetaEmbedding):

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.ordination = PCoA(self.settings)


class BatchNMDS(BatchBetaEmbedding):

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.ordination = NMDS(self.settings)


class BatchTSNE(BatchBetaEmbedding):

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.ordination = TSNE(self.settings)
