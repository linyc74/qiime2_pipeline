import pandas as pd
from typing import List
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from .tools import edit_fpath
from .template import Processor, Settings
from .embedding_core import NMDSCore, TSNECore
from .embedding_process import EmbeddingProcessTemplate


class PCA(Processor):

    NAME = 'PCA'
    XY_COLUMNS = ('PC1', 'PC2')
    LOG_PSEUDOCOUNT = True
    NORMALIZE_BY_SAMPLE_READS = False
    DATA_STRUCTURE = 'row_features'

    feature_table_tsv: str

    feature_table_df: pd.DataFrame
    sample_coordinate_df: pd.DataFrame
    proportion_explained_series: pd.Series

    def main(self, feature_table_tsv: str):
        self.feature_table_tsv = feature_table_tsv

        self.run_main_workflow()

    def run_main_workflow(self):
        self.read_feature_table_tsv()
        self.run_embedding()
        self.add_group_column()
        self.make_dstdir()
        self.write_sample_coordinate()
        self.plot_sample_coordinate()

    def read_feature_table_tsv(self):
        self.distance_matrix = pd.read_csv(
            self.distance_matrix_tsv,
            sep='\t',
            index_col=0)

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


