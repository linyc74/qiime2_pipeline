import os
import pandas as pd
import seaborn as sns
import matplotlib.axes
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import mannwhitneyu
from typing import List, Dict, Any, Tuple
from statsmodels.stats.multitest import multipletests
from .template import Processor
from .normalization import CountNormalization
from .grouping import GROUP_COLUMN, AddGroupColumn


DSTDIR_NAME = 'differential-abundance'


class DifferentialAbundance(Processor):

    taxon_table_tsv_dict: Dict[str, str]
    sample_sheet: str
    colors: List[Tuple[float, float, float, float]]
    p_value: float

    def main(
            self,
            taxon_table_tsv_dict: Dict[str, str],
            sample_sheet: str,
            colors: List[Tuple[float, float, float, float]],
            p_value: float):

        self.taxon_table_tsv_dict = taxon_table_tsv_dict
        self.sample_sheet = sample_sheet
        self.colors = colors
        self.p_value = p_value

        for taxon_level, taxon_tsv in self.taxon_table_tsv_dict.items():
            OneTaxonLevelDifferentialAbundance(self.settings).main(
                taxon_level=taxon_level,
                taxon_tsv=taxon_tsv,
                sample_sheet=self.sample_sheet,
                colors=self.colors,
                p_value=self.p_value)

        self.zip_dstdir()

    def zip_dstdir(self):
        self.call(f'tar -C "{self.outdir}" -czf "{self.outdir}/{DSTDIR_NAME}.tar.gz" {DSTDIR_NAME}')
        self.call(f'rm -r "{self.outdir}/{DSTDIR_NAME}"')


class OneTaxonLevelDifferentialAbundance(Processor):

    taxon_level: str
    taxon_tsv: str
    sample_sheet: str
    colors: List[Tuple[float, float, float, float]]
    p_value: float

    taxon_df: pd.DataFrame

    def main(
            self,
            taxon_level: str,
            taxon_tsv: str,
            sample_sheet: str,
            colors: List[Tuple[float, float, float, float]],
            p_value: float):

        self.taxon_level = taxon_level
        self.taxon_tsv = taxon_tsv
        self.sample_sheet = sample_sheet
        self.colors = colors
        self.p_value = p_value

        self.logger.info(f'Processing "{self.taxon_tsv}" at {self.taxon_level} level')
        self.read_taxon_tsv()
        self.prepare_taxon_df()
        self.mannwhitneyu_tests_and_boxplots()

    def read_taxon_tsv(self):
        self.taxon_df = pd.read_csv(self.taxon_tsv, sep='\t', index_col=0)

    def prepare_taxon_df(self):
        self.taxon_df = PrepareTaxonDf(self.settings).main(
            df=self.taxon_df,
            sample_sheet=self.sample_sheet)

    def mannwhitneyu_tests_and_boxplots(self):
        MannwhitneyuTestsAndBoxplots(self.settings).main(
            taxon_level=self.taxon_level,
            taxon_df=self.taxon_df,
            colors=self.colors,
            p_value=self.p_value)


class PrepareTaxonDf(Processor):

    df: pd.DataFrame
    sample_sheet: str

    def main(
            self,
            df: pd.DataFrame,
            sample_sheet: str) -> pd.DataFrame:

        self.df = df.copy()
        self.sample_sheet = sample_sheet

        self.normalize_counts()
        self.transpose_taxon_df()
        self.shorten_taxon_columns()
        self.add_suffix_to_duplicated_taxon_columns()
        self.add_group_column()

        return self.df

    def normalize_counts(self):
        self.df = CountNormalization(self.settings).main(
            df=self.df,
            log_pseudocount=False,
            by_sample_reads=True,
            sample_reads_unit=100)  # 100 for percentage

    def transpose_taxon_df(self):
        self.df = self.df.transpose()

    def shorten_taxon_columns(self):

        def shorten(s: str) -> str:
            """
            Seen examples where '/' is used as a separator, which causes error in the output file path
            """
            return s.replace('/', '|').split('|')[-1]

        self.df.columns = pd.Series(self.df.columns).apply(shorten)

    def add_suffix_to_duplicated_taxon_columns(self):
        self.df = AddSuffixToDuplicatedColumns(self.settings).main(self.df)

    def add_group_column(self):
        self.df = AddGroupColumn(self.settings).main(
            df=self.df,
            sample_sheet=self.sample_sheet)
        self.df[GROUP_COLUMN] = self.df[GROUP_COLUMN].astype('str')  # convert to str to be safe, int group names might cause issues


class AddSuffixToDuplicatedColumns(Processor):

    df: pd.DataFrame

    duplicated_columns: List[str]

    outdf: pd.DataFrame

    def main(self, df: pd.DataFrame) -> pd.DataFrame:
        self.df = df.copy()

        self.set_duplicated_columns()
        self.rename_duplicated_columns_and_set_outdf()

        return self.outdf

    def set_duplicated_columns(self):
        column_to_count = {}
        for c in self.df.columns:
            column_to_count.setdefault(c, 0)
            column_to_count[c] += 1

        self.duplicated_columns = [
            column for column, count in column_to_count.items()
            if count > 1
        ]

    def rename_duplicated_columns_and_set_outdf(self):
        self.outdf = pd.DataFrame(index=self.df.index)

        cumulative_count = {c: 0 for c in self.duplicated_columns}

        for column, series in self.df.items():
            if column in self.duplicated_columns:
                cumulative_count[column] += 1
                series.name = f'{column}_{cumulative_count[column]}'

            self.outdf = pd.concat([self.outdf, series], axis=1)


class MannwhitneyuTestsAndBoxplots(Processor):

    taxon_level: str
    taxon_df: pd.DataFrame
    colors: List[Tuple[float, float, float, float]]
    p_value: float

    groups: List[str]

    def main(
            self,
            taxon_level: str,
            taxon_df: pd.DataFrame,
            colors: List[Tuple[float, float, float, float]],
            p_value: float):

        self.taxon_level = taxon_level
        self.taxon_df = taxon_df
        self.colors = colors
        self.p_value = p_value

        self.plot_all()

        self.groups = self.taxon_df[GROUP_COLUMN].unique().tolist()

        for group_1, group_2 in combinations(self.groups, 2):
            self.process_group_pair(group_1=group_1, group_2=group_2)

    def plot_all(self):
        for taxon in self.taxon_df.columns:
            if taxon == GROUP_COLUMN:
                continue
            dstdir = f'{self.outdir}/{DSTDIR_NAME}/{self.taxon_level}/all'
            os.makedirs(dstdir, exist_ok=True)
            Boxplot(self.settings).main(
                data=self.taxon_df,
                x=GROUP_COLUMN,
                y=taxon,
                colors=self.colors,
                title=taxon,
                png=f'{dstdir}/{taxon}.png'
            )

    def process_group_pair(self, group_1: str, group_2: str):
        dstdir = f'{self.outdir}/{DSTDIR_NAME}/{self.taxon_level}/{group_1}-{group_2}'
        os.makedirs(dstdir, exist_ok=True)

        color_1 = self.colors[self.groups.index(group_1)]
        color_2 = self.colors[self.groups.index(group_2)]
        stats_data = []

        for taxon in self.taxon_df.columns:

            if taxon == GROUP_COLUMN:
                continue

            is_group_1 = self.taxon_df[GROUP_COLUMN] == group_1
            is_group_2 = self.taxon_df[GROUP_COLUMN] == group_2

            statistic, pvalue = mannwhitneyu(
                x=self.taxon_df.loc[is_group_1, taxon],
                y=self.taxon_df.loc[is_group_2, taxon]
            )

            if pvalue <= self.p_value:
                Boxplot(self.settings).main(
                    data=self.taxon_df[is_group_1 | is_group_2],
                    x=GROUP_COLUMN,
                    y=taxon,
                    colors=[color_1, color_2],
                    title=f'{taxon}\n$p = {pvalue:.4f}$',
                    png=f'{dstdir}/{pvalue:.4f}_{taxon}.png'
                )

            stats_data.append({
                'Taxon': taxon,
                'Mean 1 (%)': self.taxon_df.loc[is_group_1, taxon].mean(),
                'Mean 2 (%)': self.taxon_df.loc[is_group_2, taxon].mean(),
                'Statistics': statistic,
                'P value': pvalue,
            })

        self.__save_stats_data(stats_data=stats_data, dstdir=dstdir)

    def __save_stats_data(self, stats_data: List[Dict[str, Any]], dstdir: str):
        stats_df = pd.DataFrame(stats_data).sort_values(
            by='P value',
            ascending=True
        )
        rejected, pvals_corrected, _, _ = multipletests(
            stats_df['P value'],
            alpha=0.1,
            method='fdr_bh',  # Benjamini-Hochberg
            is_sorted=False,
            returnsorted=False)
        stats_df['Benjamini-Hochberg adjusted P value'] = pvals_corrected
        stats_df.to_csv(
            f'{dstdir}/Mann-Whitney-U.csv',
            index=False
        )


class Boxplot(Processor):

    WIDTH_PADDING = 1.5 / 2.54
    WIDTH_PER_GROUP = 1.25 / 2.54
    HEIGHT = 5 / 2.54
    DPI = 600
    FONT_SIZE = 6
    BOX_WIDTH = 0.5
    XLABEL = None
    YLABEL = 'Relative abundance (%)'
    LINEWIDTH = 0.5
    BOX_lINEWIDTH = 0.5
    MARKER_LINEWIDTH = 0.25
    YLIM = None

    data: pd.DataFrame
    x: str
    y: str
    colors: List[Tuple[float, float, float, float]]
    title: str
    png: str

    ax: matplotlib.axes.Axes

    def main(
            self,
            data: pd.DataFrame,
            x: str,
            y: str,
            colors: List[Tuple[float, float, float, float]],
            title: str,
            png: str):

        self.data = data
        self.x = x
        self.y = y
        self.colors = colors
        self.title = title
        self.png = png

        self.init()
        self.plot()
        self.config()
        self.save()

    def init(self):
        plt.rcParams['font.size'] = self.FONT_SIZE
        plt.rcParams['axes.linewidth'] = self.LINEWIDTH

        groups = len(self.data[self.x].unique())
        figsize = (groups * self.WIDTH_PER_GROUP + self.WIDTH_PADDING, self.HEIGHT)

        plt.figure(figsize=figsize)

    def plot(self):
        self.ax = sns.boxplot(
            data=self.data,
            x=self.x,
            y=self.y,
            hue=self.x,
            palette=self.colors,
            width=self.BOX_WIDTH,
            linewidth=self.BOX_lINEWIDTH,
            dodge=False,  # to align the boxes on the x axis
        )
        self.ax = sns.stripplot(
            data=self.data,
            x=self.x,
            y=self.y,
            hue=self.x,
            palette=self.colors,
            linewidth=self.MARKER_LINEWIDTH,
        )

    def config(self):
        self.ax.set_title(self.title)
        self.ax.set(xlabel=self.XLABEL, ylabel=self.YLABEL)
        plt.gca().xaxis.set_tick_params(width=self.LINEWIDTH)
        plt.gca().yaxis.set_tick_params(width=self.LINEWIDTH)
        plt.ylim(self.YLIM)
        plt.legend().remove()

    def save(self):
        plt.tight_layout()
        plt.savefig(self.png, dpi=self.DPI)
        plt.close()
