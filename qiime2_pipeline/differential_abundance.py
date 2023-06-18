import os
import pandas as pd
import seaborn as sns
import matplotlib.axes
import matplotlib.pyplot as plt
from itertools import combinations
from typing import List, Tuple, Dict
from scipy.stats import mannwhitneyu
from .template import Processor
from .grouping import AddGroupColumn
from .normalization import CountNormalization


GROUP_COLUMN = AddGroupColumn.GROUP_COLUMN
DSTDIR_NAME = 'differential-abundance'
FIGSIZE = (4 / 2.54, 5 / 2.54)
DPI = 600
FONT_SIZE = 7
BOX_WIDTH = 0.5
PALETTE = 'CMRmap'
XLABEL = None
YLABEL = 'Relative abundance (%)'
MARKER_LINEWIDTH = 0.5
YLIM = None


class DifferentialAbundance(Processor):

    taxon_table_tsv_dict: Dict[str, str]
    group_keywords: List[str]

    def main(
            self,
            taxon_table_tsv_dict: Dict[str, str],
            group_keywords: List[str]):

        self.taxon_table_tsv_dict = taxon_table_tsv_dict
        self.group_keywords = group_keywords

        for taxon_level, taxon_tsv in self.taxon_table_tsv_dict.items():
            OneTaxonLevelDifferentialAbundance(self.settings).main(
                taxon_level=taxon_level,
                taxon_tsv=taxon_tsv,
                group_keywords=self.group_keywords)

        self.zip_dstdir()

    def zip_dstdir(self):
        dirpath = f'{self.outdir}/{DSTDIR_NAME}'
        self.call(f'zip -r {dirpath}.zip {dirpath} >> {self.outdir}/zip.log')
        self.call(f'rm -r {dirpath}')


class OneTaxonLevelDifferentialAbundance(Processor):

    taxon_level: str
    taxon_tsv: str
    group_keywords: List[str]

    taxon_df: pd.DataFrame

    def main(
            self,
            taxon_level: str,
            taxon_tsv: str,
            group_keywords: List[str]):

        self.taxon_level = taxon_level
        self.taxon_tsv = taxon_tsv
        self.group_keywords = group_keywords

        self.read_taxon_tsv()
        self.prepare_taxon_df()
        self.mannwhitneyu_tests_and_boxplots()

    def read_taxon_tsv(self):
        self.taxon_df = pd.read_csv(self.taxon_tsv, sep='\t', index_col=0)

    def prepare_taxon_df(self):
        self.taxon_df = PrepareTaxonDf(self.settings).main(
            df=self.taxon_df,
            group_keywords=self.group_keywords)

    def mannwhitneyu_tests_and_boxplots(self):
        for group_1, group_2 in combinations(self.group_keywords, 2):

            dstdir = f'{self.outdir}/{DSTDIR_NAME}/{self.taxon_level}/{group_1}-{group_2}'
            os.makedirs(dstdir, exist_ok=True)
            stats_data = []

            for taxon in self.taxon_df.columns:

                if taxon == GROUP_COLUMN:
                    continue

                statistic, pvalue = MannwhitneyuTestAndBoxplot(self.settings).main(
                    df=self.taxon_df,
                    group_1=group_1,
                    group_2=group_2,
                    taxon=taxon,
                    dstdir=dstdir)

                stats_data.append({
                    'Taxon': taxon,
                    'Statistics': statistic,
                    'P value': pvalue
                })

            pd.DataFrame(stats_data).sort_values(
                by='P value',
                ascending=True
            ).to_csv(
                f'{dstdir}/Mann-Whitney-U.csv',
                index=False
            )


class PrepareTaxonDf(Processor):

    df: pd.DataFrame
    group_keywords: List[str]

    def main(
            self,
            df: pd.DataFrame,
            group_keywords: List[str]) -> pd.DataFrame:

        self.df = df.copy()
        self.group_keywords = group_keywords

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
            return s.split('|')[-1]
        self.df.columns = pd.Series(self.df.columns).apply(shorten)

    def add_suffix_to_duplicated_taxon_columns(self):
        self.df = AddSuffixToDuplicatedColumns(self.settings).main(self.df)

    def add_group_column(self):
        self.df = AddGroupColumn(self.settings).main(
            df=self.df,
            group_keywords=self.group_keywords)


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

        for column, series in self.df.iteritems():
            if column in self.duplicated_columns:
                cumulative_count[column] += 1
                series.name = f'{column}_{cumulative_count[column]}'

            self.outdf = pd.concat([self.outdf, series], axis=1)


class MannwhitneyuTestAndBoxplot(Processor):

    df: pd.DataFrame
    group_1: str
    group_2: str
    taxon: str
    dstdir: str

    statistic: float
    pvalue: float

    def main(
            self,
            df: pd.DataFrame,
            group_1: str,
            group_2: str,
            taxon: str,
            dstdir: str) -> Tuple[float, float]:

        self.df = df
        self.group_1 = group_1
        self.group_2 = group_2
        self.taxon = taxon
        self.dstdir = dstdir

        self.mannwhitneyu_test()
        self.boxplot()

        return self.statistic, self.pvalue

    def mannwhitneyu_test(self):
        is_group_1 = self.df[GROUP_COLUMN] == self.group_1
        is_group_2 = self.df[GROUP_COLUMN] == self.group_2

        res = mannwhitneyu(
            x=self.df.loc[is_group_1, self.taxon],
            y=self.df.loc[is_group_2, self.taxon]
        )

        self.statistic, self.pvalue = res.statistic, res.pvalue

    def boxplot(self):
        isin_groups = self.df[GROUP_COLUMN].isin([self.group_1, self.group_2])
        Boxplot(self.settings).main(
            data=self.df[isin_groups],
            x=GROUP_COLUMN,
            y=self.taxon,
            title=f'p = {self.pvalue:.4f}',
            dstdir=self.dstdir)


class Boxplot(Processor):

    data: pd.DataFrame
    x: str
    y: str
    title: str
    dstdir: str

    ax: matplotlib.axes.Axes

    def main(self, data: pd.DataFrame, x: str, y: str, title: str, dstdir: str):
        self.data = data
        self.x = x
        self.y = y
        self.title = title
        self.dstdir = dstdir

        self.init()
        self.plot()
        self.config()
        self.save()

    def init(self):
        plt.figure(figsize=FIGSIZE)
        plt.rc('font', size=FONT_SIZE)

    def plot(self):
        self.ax = sns.boxplot(
            data=self.data,
            x=self.x,
            y=self.y,
            palette=PALETTE,
            width=BOX_WIDTH
        )
        self.ax = sns.stripplot(
            data=self.data,
            x=self.x,
            y=self.y,
            palette=PALETTE,
            linewidth=MARKER_LINEWIDTH
        )

    def config(self):
        self.ax.set_title(self.title)
        self.ax.set(xlabel=XLABEL, ylabel=YLABEL)
        plt.ylim(YLIM)

    def save(self):
        plt.tight_layout()
        plt.savefig(f'{self.dstdir}/{self.y}.png', dpi=DPI)
        plt.close()
