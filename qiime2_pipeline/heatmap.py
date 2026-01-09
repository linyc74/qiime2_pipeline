import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple, List
from .utils import edit_fpath
from .template import Processor
from .normalization import CountNormalization
from .grouping import TagGroupNamesOnSampleColumns


class PlotHeatmaps(Processor):

    DSTDIR_NAME = 'heatmap'

    tsvs: List[str]
    heatmap_read_fraction: float
    sample_sheet: str

    dstdir: str

    def main(
            self,
            tsvs: List[str],
            heatmap_read_fraction: float,
            sample_sheet: str):

        self.tsvs = tsvs
        self.heatmap_read_fraction = heatmap_read_fraction
        self.sample_sheet = sample_sheet

        self.make_dstdir()
        for tsv in self.tsvs:
            self.plot_one_heatmap(tsv=tsv)

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def plot_one_heatmap(self, tsv: str):
        PlotOneHeatmap(self.settings).main(
            tsv=tsv,
            heatmap_read_fraction=self.heatmap_read_fraction,
            sample_sheet=self.sample_sheet,
            dstdir=self.dstdir)


class PlotOneHeatmap(Processor):

    LOG_PSEUDOCOUNT = True
    NORMALIZE_BY_SAMPLE_READS = True

    tsv: str
    heatmap_read_fraction: float
    sample_sheet: str
    dstdir: str

    df: pd.DataFrame

    def main(
            self,
            tsv: str,
            heatmap_read_fraction: float,
            sample_sheet: str,
            dstdir: str):

        self.tsv = tsv
        self.heatmap_read_fraction = heatmap_read_fraction
        self.sample_sheet = sample_sheet
        self.dstdir = dstdir

        self.read_tsv()
        self.filter_by_cumulative_reads()
        self.count_normalization()
        self.clustermap()

    def read_tsv(self):
        self.df = pd.read_csv(self.tsv, sep='\t', index_col=0)

    def filter_by_cumulative_reads(self):
        self.df = FilterByCumulativeReads(self.settings).main(
            df=self.df,
            heatmap_read_fraction=self.heatmap_read_fraction)

    def count_normalization(self):
        self.df = CountNormalization(self.settings).main(
            df=self.df,
            log_pseudocount=self.LOG_PSEUDOCOUNT,
            by_sample_reads=self.NORMALIZE_BY_SAMPLE_READS)

    def clustermap(self):
        for cluster_columns, suffix in [
            (True, '-sample-clustered'),
            (False, '-sample-unclustered')
        ]:
            output_prefix = edit_fpath(
                fpath=self.tsv,
                old_suffix='.tsv',
                new_suffix=suffix,
                dstdir=self.dstdir)

            Clustermap(self.settings).main(
                data=self.df,
                sample_sheet=self.sample_sheet,
                cluster_columns=cluster_columns,
                output_prefix=output_prefix)


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

    COLORMAP = 'PuBu'
    Y_LABEL_CHAR_WIDTH = 0.14 / 2.54
    X_LABEL_CHAR_WIDTH = 0.14 / 2.54
    CELL_WIDTH = 0.4 / 2.54
    CELL_HEIGHT = 0.4 / 2.54
    DENDROGRAM_SIZE = 1.0 / 2.54
    COLORBAR_WIDTH = 0.01
    COLORBAR_HORIZONTAL_POSITION = 1.
    FONTSIZE = 7
    LINE_WIDTH = 0.5
    DPI = 600

    data: pd.DataFrame
    sample_sheet: str
    cluster_columns: bool
    output_prefix: str

    x_label_padding: float
    y_label_padding: float
    figsize: Tuple[float, float]
    grid: sns.matrix.ClusterGrid

    def main(
            self,
            data: pd.DataFrame,
            sample_sheet: str,
            cluster_columns: bool,
            output_prefix: str):

        self.data = data.copy()
        self.sample_sheet = sample_sheet
        self.cluster_columns = cluster_columns
        self.output_prefix = output_prefix

        self.tag_group_names_on_sample_columns()
        self.shorten_taxon_names_for_publication()
        self.set_figsize()
        self.clustermap()
        self.config_clustermap()
        self.save_fig()
        self.save_csv()

    def tag_group_names_on_sample_columns(self):
        self.data = TagGroupNamesOnSampleColumns(self.settings).main(
            df=self.data,
            sample_sheet=self.sample_sheet)

    def shorten_taxon_names_for_publication(self):
        if self.settings.for_publication:
            self.data.index = self.data.index.str.split('|').str[-1]

    def set_figsize(self):
        self.__set_x_y_label_padding()
        w = (len(self.data.columns) * self.CELL_WIDTH) + self.y_label_padding
        h = (len(self.data.index) * self.CELL_HEIGHT) + self.x_label_padding
        self.figsize = (w, h)

    def __set_x_y_label_padding(self):
        max_x_label_length = pd.Series(self.data.columns).apply(len).max()
        self.x_label_padding = max_x_label_length * self.X_LABEL_CHAR_WIDTH + self.DENDROGRAM_SIZE

        max_y_label_length = pd.Series(self.data.index).apply(len).max()
        self.y_label_padding = max_y_label_length * self.Y_LABEL_CHAR_WIDTH + self.DENDROGRAM_SIZE

    def clustermap(self):
        plt.rcParams['font.size'] = self.FONTSIZE
        plt.rcParams['axes.linewidth'] = self.LINE_WIDTH

        row_cluster = True if len(self.data) > 1 else False  # do not cluster if only one row

        w, h = self.figsize
        dendrogram_ratio = (self.DENDROGRAM_SIZE / w, self.DENDROGRAM_SIZE / h)

        self.grid = sns.clustermap(
            data=self.data,
            cmap=self.COLORMAP,
            figsize=self.figsize,
            xticklabels=True,  # include every x label
            yticklabels=True,  # include every y label
            row_cluster=row_cluster,
            col_cluster=self.cluster_columns,
            dendrogram_ratio=dendrogram_ratio,
            linewidth=0.25)
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
            labelrotation=0,
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

    def save_fig(self):
        # must use grid.savefig(), but not plt.savefig()
        # plt.savefig() crops out the colorbar
        plt.tight_layout()
        dpi = self.__downsize_dpi_if_too_large()
        for ext in ['pdf', 'png']:
            self.grid.savefig(f'{self.output_prefix}.{ext}', dpi=dpi)
        plt.close()

    def __downsize_dpi_if_too_large(self) -> int:
        longer_side = max(self.figsize)
        dpi = self.DPI
        while longer_side * dpi >= 2**15:  # 2^16 is the limit of matplotlib, but 2^15 is safer after some tests
            dpi = int(dpi/2)  # downsize
        return dpi

    def save_csv(self):
        self.data.to_csv(f'{self.output_prefix}.tsv', sep='\t', index=True)
