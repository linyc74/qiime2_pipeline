import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Tuple
from abc import ABC, abstractmethod
from .tools import edit_fpath
from .template import Processor
from .grouping import AddGroupColumn
from matplotlib.axes import Axes


class EmbeddingProcessTemplate(Processor, ABC):

    NAME: str
    XY_COLUMNS: Tuple[str, str]
    DSTDIR_NAME: str
    GROUP_COLUMN: str = AddGroupColumn.GROUP_COLUMN

    tsv: str
    group_keywords: List[str]

    df: pd.DataFrame
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

        self.preprocessing()
        self.embedding()

        self.add_group_column()
        self.make_dstdir()
        self.write_sample_coordinate()
        self.plot_sample_coordinate()

    def read_tsv(self):
        self.df = pd.read_csv(
            self.tsv,
            sep='\t',
            index_col=0)

    @abstractmethod
    def preprocessing(self):
        pass

    @abstractmethod
    def embedding(self):
        pass

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def add_group_column(self):
        self.sample_coordinate_df = AddGroupColumn(self.settings).main(
            df=self.sample_coordinate_df,
            group_keywords=self.group_keywords)

    def write_sample_coordinate(self):
        tsv = self.__get_sample_coordinate_fpath(suffix='.tsv')
        self.sample_coordinate_df.to_csv(tsv, sep='\t')

    def plot_sample_coordinate(self):
        output_prefix = self.__get_sample_coordinate_fpath(suffix='')
        ScatterPlot(self.settings).main(
            sample_coordinate_df=self.sample_coordinate_df,
            x_column=self.XY_COLUMNS[0],
            y_column=self.XY_COLUMNS[1],
            hue_column=self.GROUP_COLUMN,
            output_prefix=output_prefix)

    def __get_sample_coordinate_fpath(self, suffix: str) -> str:
        name = self.NAME.lower().replace('-', '')
        return edit_fpath(
            fpath=self.tsv,
            old_suffix='.tsv',
            new_suffix=f'-{name}-sample-coordinate{suffix}',
            dstdir=self.dstdir
        )


class ScatterPlot(Processor):

    FIGSIZE = (8, 8)
    DPI = 600

    sample_coordinate_df: pd.DataFrame
    x_column: str
    y_column: str
    group_column: str
    output_prefix: str

    ax: Axes

    def main(
            self,
            sample_coordinate_df: pd.DataFrame,
            x_column: str,
            y_column: str,
            hue_column: str,
            output_prefix: str):

        self.sample_coordinate_df = sample_coordinate_df
        self.x_column = x_column
        self.y_column = y_column
        self.group_column = hue_column
        self.output_prefix = output_prefix

        self.init_figure()
        self.scatterplot()
        self.label_points()
        self.save_figure()

    def init_figure(self):
        plt.figure(figsize=self.FIGSIZE, dpi=self.DPI)

    def scatterplot(self):
        self.ax = sns.scatterplot(
            data=self.sample_coordinate_df,
            x=self.x_column,
            y=self.y_column,
            hue=self.group_column)

    def label_points(self):
        df = self.sample_coordinate_df
        for sample_name in df.index:
            self.ax.text(
                x=df.loc[sample_name, self.x_column],
                y=df.loc[sample_name, self.y_column],
                s=sample_name
            )

    def save_figure(self):
        for ext in ['pdf', 'png']:
            plt.savefig(f'{self.output_prefix}.{ext}', dpi=self.DPI)
        plt.close()
