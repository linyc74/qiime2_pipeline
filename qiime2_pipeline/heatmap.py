import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple, List
from .tools import edit_fpath
from .template import Processor


class PlotHeatmaps(Processor):

    DSTDIR_NAME = 'heatmap'

    tsvs: List[str]
    heatmap_read_fraction: float

    dstdir: str

    def main(
            self,
            tsvs: List[str],
            heatmap_read_fraction: float):

        self.tsvs = tsvs
        self.heatmap_read_fraction = heatmap_read_fraction

        self.make_dstdir()
        for tsv in self.tsvs:
            self.plot_one(tsv=tsv)

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def plot_one(self, tsv: str):
        data = pd.read_csv(tsv, sep='\t', index_col=0)
        data = FilterByCumulativeReads(self.settings).main(
            df=data,
            heatmap_read_fraction=self.heatmap_read_fraction
        )
        output_prefix = edit_fpath(
            fpath=tsv,
            old_suffix='.tsv',
            new_suffix='',
            dstdir=self.dstdir
        )
        Clustermap(self.settings).main(
            data=data,
            output_prefix=output_prefix
        )


class FilterByCumulativeReads(Processor):

    ROW_SUM = 'Row Sum'
    CUMULATIVE_SUM = 'Cumulative Sum'

    df: pd.DataFrame
    heatmap_read_fraction: float

    def main(
            self,
            df: pd.DataFrame,
            heatmap_read_fraction: float) -> pd.DataFrame:

        self.df = df
        self.heatmap_read_fraction = heatmap_read_fraction

        self.sum_each_row()
        self.sort_by_sum()
        self.cumulate_sum()
        self.divide_by_total()
        self.filter()
        self.clean_up()

        return self.df

    def sum_each_row(self):
        self.df[self.ROW_SUM] = np.sum(self.df.to_numpy(), axis=1)

    def sort_by_sum(self):
        self.df = self.df.sort_values(
            by=self.ROW_SUM,
            ascending=False)

    def cumulate_sum(self):
        l = []
        c = 0
        for s in self.df[self.ROW_SUM]:
            c += s
            l.append(c)
        self.df[self.CUMULATIVE_SUM] = l

    def divide_by_total(self):
        total = self.df[self.CUMULATIVE_SUM][-1]
        self.df[self.CUMULATIVE_SUM] = self.df[self.CUMULATIVE_SUM] / total

    def filter(self):
        assert 0 < self.heatmap_read_fraction < 1
        n_rows = len(self.df)
        for i, v in enumerate(self.df[self.CUMULATIVE_SUM]):
            if v > self.heatmap_read_fraction:
                n_rows = i
                break
        self.df = self.df.iloc[:n_rows, ]

    def clean_up(self):
        self.df = self.df.drop(
            columns=[self.ROW_SUM, self.CUMULATIVE_SUM]
        )


class Clustermap(Processor):

    DSTDIR_NAME = 'heatmap'
    CLUSTER_COLUMNS = True
    COLORMAP = 'PuBu'
    Y_LABEL_CHAR_WIDTH = 0.08
    VERTICAL_PADDING = 2
    CELL_WIDTH = 0.3
    CELL_HEIGHT = 0.3
    ROW_TREE_WIDTH = 1
    COL_TREE_HEIGHT = 1
    COLORBAR_WIDTH = 0.01
    COLORBAR_HEIGHT = 0.1
    COLORBAR_HORIZONTAL_POSITION = 1.

    data: pd.DataFrame
    output_prefix: str

    horizontal_paddng: float
    figsize: Tuple[float, float]
    dendrogram_ratio: Tuple[float, float]
    grid: sns.matrix.ClusterGrid

    def main(self, data: pd.DataFrame, output_prefix: str):
        self.data = data
        self.output_prefix = output_prefix

        self.set_figsize()
        self.set_dendrogram_ratio()
        self.clustermap()
        self.config_clustermap()
        self.save_pdf()
        self.save_csv()

    def set_figsize(self):
        self.__set_horizontal_padding()
        w = (len(self.data.columns) * self.CELL_WIDTH) + self.horizontal_paddng
        h = (len(self.data.index) * self.CELL_HEIGHT) + self.VERTICAL_PADDING
        self.figsize = (w, h)

    def __set_horizontal_padding(self):
        max_y_label_length = pd.Series(self.data.index).apply(len).max()
        self.horizontal_paddng = max_y_label_length * self.Y_LABEL_CHAR_WIDTH

    def set_dendrogram_ratio(self):
        heatmap_width = len(self.data.columns) * self.CELL_WIDTH
        heatmap_height = len(self.data.index) * self.CELL_HEIGHT
        row_tree_ratio = self.ROW_TREE_WIDTH / heatmap_width
        col_tree_ratio = self.COL_TREE_HEIGHT / heatmap_height
        self.dendrogram_ratio = (row_tree_ratio, col_tree_ratio)

    def clustermap(self):
        self.grid = sns.clustermap(
            data=self.data,
            cmap=self.COLORMAP,
            figsize=self.figsize,
            xticklabels=True,  # include every x label
            yticklabels=True,  # include every y label
            col_cluster=self.CLUSTER_COLUMNS,
            dendrogram_ratio=self.dendrogram_ratio,
            linewidth=0.5)
        self.__set_plotted_data()

    def __set_plotted_data(self):
        self.data = self.grid.__dict__['data2d']

    def config_clustermap(self):
        self.__set_x_y_axes()
        self.__set_colorbar()

    def __set_x_y_axes(self):
        heatmap = self.grid.ax_heatmap

        n_rows = len(self.data)
        heatmap.set_ylim(n_rows, 0)  # first and last row not be chopped half

        heatmap.tick_params(
            axis='x',
            bottom=True,
            top=False,
            labelbottom=True,
            labeltop=False,
            labelrotation=90
        )
        heatmap.tick_params(
            axis='y',
            left=False,
            right=True,
            labelleft=False,
            labelright=True,
            labelrotation=0
        )

    def __set_colorbar(self):
        colorbar = self.grid.cax
        p = colorbar.get_position()
        colorbar.set_position([
            self.COLORBAR_HORIZONTAL_POSITION,  # left
            p.y0,  # bottom
            self.COLORBAR_WIDTH,  # width
            p.height  # height
        ])

    def save_pdf(self):
        # must use grid.savefig(), but not plt.savefig()
        # plt.savefig() crops out the colorbar
        self.grid.savefig(f'{self.output_prefix}.pdf')
        plt.close()

    def save_csv(self):
        self.data.to_csv(f'{self.output_prefix}.tsv', sep='\t', index=True)
