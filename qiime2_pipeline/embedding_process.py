import os
import pandas as pd
from typing import List, Tuple
from abc import abstractmethod
from .tools import edit_fpath
from .template import Processor
from .grouping import AddGroupColumn
from .embedding_core import ScatterPlot


class EmbeddingProcessTemplate(Processor):

    NAME: str
    XY_COLUMNS: Tuple[str, str]
    DSTDIR_NAME: str
    GROUP_COLUMN: str = AddGroupColumn.GROUP_COLUMN

    tsv: str
    group_keywords: List[str]

    distance_matrix: pd.DataFrame
    sample_coordinate_df: pd.DataFrame
    dstdir: str

    @abstractmethod
    def main(
            self,
            tsv: str,
            group_keywords: List[str]):
        pass

    def run_main_workflow(self):
        self.read_tsv()
        self.run_embedding()
        self.add_group_column()
        self.make_dstdir()
        self.write_sample_coordinate()
        self.plot_sample_coordinate()

    def read_tsv(self):
        self.distance_matrix = pd.read_csv(
            self.tsv,
            sep='\t',
            index_col=0)

    @abstractmethod
    def run_embedding(self):
        pass

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
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
            fpath=self.tsv,
            old_suffix='.tsv',
            new_suffix=f'-{name}-sample-coordinate.{ext}',
            dstdir=self.dstdir
        )
