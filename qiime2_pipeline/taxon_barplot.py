import os
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple, Dict
from .template import Processor
from .normalization import CountNormalization


class PlotTaxonBarplots(Processor):

    DSTDIR_NAME = 'taxon-barplot'

    taxon_table_tsv_dict: Dict[str, str]
    n_taxa: int

    dstdir: str

    def main(
            self,
            taxon_table_tsv_dict: Dict[str, str],
            n_taxa: int):

        self.taxon_table_tsv_dict = taxon_table_tsv_dict
        self.n_taxa = n_taxa

        self.make_dstdir()
        for level, tsv in self.taxon_table_tsv_dict.items():
            self.plot_one_taxon_barplot(level=level, tsv=tsv)

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def plot_one_taxon_barplot(
            self,
            level: str,
            tsv: str):

        PlotOneTaxonBarplot(self.settings).main(
            taxon_level=level,
            taxon_table_tsv=tsv,
            n_taxa=self.n_taxa,
            dstdir=self.dstdir)


class PlotOneTaxonBarplot(Processor):

    taxon_level: str
    taxon_table_tsv: str
    n_taxa: int
    dstdir: str

    df: pd.DataFrame

    def main(
            self,
            taxon_level: str,
            taxon_table_tsv: str,
            n_taxa: int,
            dstdir: str):

        self.taxon_level = taxon_level
        self.taxon_table_tsv = taxon_table_tsv
        self.n_taxa = n_taxa
        self.dstdir = dstdir

        self.read_tsv()
        self.pool_minor_taxa()
        self.percentage_normalization()
        self.save_tsv()
        self.parcentage_barplot()

    def read_tsv(self):
        self.df = pd.read_csv(self.taxon_table_tsv, sep='\t', index_col=0)

    def pool_minor_taxa(self):
        self.df = PoolMinorFeatures(self.settings).main(
            df=self.df,
            n_major_features=self.n_taxa)

    def percentage_normalization(self):
        self.df = CountNormalization(self.settings).main(
            df=self.df,
            log_pseudocount=False,
            by_sample_reads=True,
            sample_reads_unit=100)

    def save_tsv(self):
        self.df.to_csv(f'{self.dstdir}/{self.taxon_level}-barplot.tsv', sep='\t')

    def parcentage_barplot(self):
        PercentageBarplot(self.settings).main(
            data=self.df,
            title=self.taxon_level,
            output_png=f'{self.dstdir}/{self.taxon_level}-barplot.png'
        )


class PoolMinorFeatures(Processor):

    POOLED_FEATURE_NAME = 'Others'
    ROW_SUM = 'Row Sum'

    df: pd.DataFrame
    n_major_features: int

    def main(
            self,
            df: pd.DataFrame,
            n_major_features: int) -> pd.DataFrame:

        self.df = df
        self.n_major_features = n_major_features

        self.sum_each_row()
        self.sort_by_sum()
        self.pool_minor_features()
        self.clean_up()

        return self.df

    def sum_each_row(self):
        self.df[self.ROW_SUM] = np.sum(self.df.to_numpy(), axis=1)

    def sort_by_sum(self):
        self.df = self.df.sort_values(
            by=self.ROW_SUM,
            ascending=False)

    def pool_minor_features(self):
        if len(self.df) <= self.n_major_features:
            return
        minor_feature_df = self.df.iloc[self.n_major_features:]  # extract minor features
        self.df = self.df.iloc[0:self.n_major_features]  # remove from main df
        self.df.loc[self.POOLED_FEATURE_NAME] = np.sum(minor_feature_df, axis=0)  # sum -> new row

    def clean_up(self):
        self.df = self.df.drop(columns=self.ROW_SUM)


class PercentageBarplot(Processor):

    DPI = 300
    Y_LABEL = 'Percentage'
    X_LABEL = 'Sample'
    X_LABEL_CHAR_WIDTH = 0.08
    LEGEND_CHAR_WIDTH = 0.08
    BAR_WIDTH = 0.3
    FEATURE_HEIGHT = 0.3
    COLORS = [
        'lavenderblush',
        'midnightblue',
        'royalblue',
        'cornflowerblue',
        'purple',
        'palevioletred',
        'mediumvioletred',
        'moccasin',
        'firebrick',
        'silver',
        'rebeccapurple',
        'turquoise',
        'yellow',
        'crimson',
        'orangered',
        'darkgreen',
    ]

    data: pd.DataFrame
    title: str
    output_png: str

    figsize: Tuple[float, float]
    figure: plt.Figure

    def main(
            self,
            data: pd.DataFrame,
            title: str,
            output_png: str):

        self.data = data
        self.title = title
        self.output_png = output_png

        self.set_figsize()
        self.init_figure()
        self.plot()
        self.config_and_save()

    def set_figsize(self):
        self.__set_paddings()
        w = (len(self.data.columns) * self.BAR_WIDTH) + self.horizontal_padding
        h = (len(self.data.index) * self.FEATURE_HEIGHT) + self.vertical_padding
        self.figsize = (w, h)

    def __set_paddings(self):
        max_x_label_length = pd.Series(self.data.columns).apply(len).max()
        self.vertical_padding = max_x_label_length * self.X_LABEL_CHAR_WIDTH

        max_legend_length = pd.Series(self.data.index).apply(len).max()
        self.horizontal_padding = max_legend_length * self.LEGEND_CHAR_WIDTH

    def init_figure(self):
        self.figure = plt.figure(figsize=self.figsize, dpi=self.DPI)

    def plot(self):
        random.seed(1)  # reproducible color
        bottom = np.zeros(shape=len(self.data.columns))
        for i, feature in enumerate(self.data.index):
            color = self.COLORS[i] \
                if i < len(self.COLORS) \
                else get_random_color()
            plt.bar(
                x=self.data.columns,
                height=self.data.loc[feature],
                bottom=bottom,
                color=color,
                width=0.85,
                edgecolor='black'
            )
            bottom += self.data.loc[feature]

    def config_and_save(self):
        plt.xlim(left=-1, right=len(self.data.columns))
        plt.xticks(rotation=90)
        plt.title(self.title)
        plt.xlabel(self.X_LABEL)
        plt.ylabel(self.Y_LABEL)
        plt.legend(self.data.index, bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.savefig(self.output_png)


def get_random_color() -> Tuple[float]:
    return tuple(random.randrange(150, 256) / 256 for _ in range(3))
