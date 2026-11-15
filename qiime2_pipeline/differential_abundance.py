import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.axes
import matplotlib.figure
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import mannwhitneyu, wilcoxon
from typing import List, Dict, Any, Tuple, Optional
from statsmodels.stats.multitest import multipletests
from .template import Processor
from .normalization import CountNormalization
from .grouping import GROUP_COLUMN, AddGroupColumn


DSTDIR_NAME = 'differential-abundance'


class DifferentialAbundance(Processor):

    taxon_table_tsv_dict: Dict[str, str]
    sample_sheet: str
    subject_column: Optional[str]
    colors: List[Tuple[float, float, float, float]]
    p_value: float

    def main(
            self,
            taxon_table_tsv_dict: Dict[str, str],
            sample_sheet: str,
            subject_column: Optional[str],
            colors: List[Tuple[float, float, float, float]],
            p_value: float,
            min_abundance_per_group: float):

        self.taxon_table_tsv_dict = taxon_table_tsv_dict
        self.sample_sheet = sample_sheet
        self.subject_column = subject_column
        self.colors = colors
        self.p_value = p_value
        self.min_abundance_per_group = min_abundance_per_group

        for taxon_level, taxon_tsv in self.taxon_table_tsv_dict.items():
            OneTaxonLevelDifferentialAbundance(self.settings).main(
                taxon_level=taxon_level,
                taxon_tsv=taxon_tsv,
                sample_sheet=self.sample_sheet,
                subject_column=self.subject_column,
                colors=self.colors,
                p_value=self.p_value,
                min_abundance_per_group=self.min_abundance_per_group)

        self.zip_dstdir()

    def zip_dstdir(self):
        self.call(f'tar -C "{self.outdir}" -czf "{self.outdir}/{DSTDIR_NAME}.tar.gz" {DSTDIR_NAME}')
        self.call(f'rm -r "{self.outdir}/{DSTDIR_NAME}"')


class OneTaxonLevelDifferentialAbundance(Processor):

    taxon_level: str
    taxon_tsv: str
    sample_sheet: str
    subject_column: Optional[str]
    colors: List[Tuple[float, float, float, float]]
    p_value: float

    taxon_df: pd.DataFrame

    def main(
            self,
            taxon_level: str,
            taxon_tsv: str,
            sample_sheet: str,
            subject_column: Optional[str],
            colors: List[Tuple[float, float, float, float]],
            p_value: float,
            min_abundance_per_group: float):

        self.taxon_level = taxon_level
        self.taxon_tsv = taxon_tsv
        self.sample_sheet = sample_sheet
        self.subject_column = subject_column
        self.colors = colors
        self.p_value = p_value
        self.min_abundance_per_group = min_abundance_per_group

        self.logger.info(f'Processing "{self.taxon_tsv}" at {self.taxon_level} level')
        self.read_taxon_tsv()
        self.prepare_taxon_df()
        self.statistics_and_boxplots()

    def read_taxon_tsv(self):
        self.taxon_df = pd.read_csv(self.taxon_tsv, sep='\t', index_col=0)

    def prepare_taxon_df(self):
        self.taxon_df = PrepareTaxonDf(self.settings).main(
            df=self.taxon_df,
            sample_sheet=self.sample_sheet,
            subject_column=self.subject_column)

    def statistics_and_boxplots(self):
        StatisticsAndBoxplots(self.settings).main(
            taxon_level=self.taxon_level,
            taxon_df=self.taxon_df,
            colors=self.colors,
            p_value=self.p_value,
            min_abundance_per_group=self.min_abundance_per_group,
            subject_column=self.subject_column)


class PrepareTaxonDf(Processor):

    df: pd.DataFrame
    sample_sheet: str
    subject_column: Optional[str]

    def main(
            self,
            df: pd.DataFrame,
            sample_sheet: str,
            subject_column: Optional[str]) -> pd.DataFrame:

        self.df = df.copy()
        self.sample_sheet = sample_sheet
        self.subject_column = subject_column

        self.normalize_counts()
        self.transpose_taxon_df()
        self.shorten_taxon_columns()
        self.add_suffix_to_duplicated_taxon_columns()
        self.add_group_column()
        if self.subject_column is not None:
            self.add_subject_column()

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

    def add_subject_column(self):
        sample_sheet_df = pd.read_csv(self.sample_sheet, index_col=0)[self.subject_column]  # only take the subject column
        self.df = self.df.merge(
            right=sample_sheet_df,
            how='left',
            left_index=True,
            right_index=True
        )
        cols = list(self.df.columns)
        cols = [cols[0]] + [cols[-1]] + cols[1:-1]  # move the subject column to the second position
        self.df = self.df[cols]


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


class StatisticsAndBoxplots(Processor):

    taxon_level: str
    taxon_df: pd.DataFrame
    colors: List[Tuple[float, float, float, float]]
    p_value: float
    min_abundance_per_group: float
    subject_column: Optional[str]

    groups: List[str]

    def main(
            self,
            taxon_level: str,
            taxon_df: pd.DataFrame,
            colors: List[Tuple[float, float, float, float]],
            p_value: float,
            min_abundance_per_group: float,
            subject_column: Optional[str]):

        self.taxon_level = taxon_level
        self.taxon_df = taxon_df
        self.colors = colors
        self.p_value = p_value
        self.min_abundance_per_group = min_abundance_per_group
        self.subject_column = subject_column

        self.plot_all_taxa()

        self.groups = self.taxon_df[GROUP_COLUMN].unique().tolist()

        for group1, group2 in combinations(self.groups, 2):
            self.process_group_pair(group1=group1, group2=group2)

    def plot_all_taxa(self):
        for taxon in self.taxon_df.columns:
            if taxon in [GROUP_COLUMN, self.subject_column]:
                continue
            dstdir = f'{self.outdir}/{DSTDIR_NAME}/{self.taxon_level}/all'
            os.makedirs(dstdir, exist_ok=True)
            Boxplot(self.settings).main(
                data=self.taxon_df,
                x=GROUP_COLUMN,
                y=taxon,
                units=self.subject_column,
                colors=self.colors,
                title=taxon,
                png=f'{dstdir}/{taxon}.png'
            )

    def process_group_pair(self, group1: str, group2: str):
        dstdir = f'{self.outdir}/{DSTDIR_NAME}/{self.taxon_level}/{group1}-{group2}'
        os.makedirs(dstdir, exist_ok=True)

        stats_data = []

        for taxon in self.taxon_df.columns:

            if taxon in [GROUP_COLUMN, self.subject_column]:
                continue

            if self.subject_column is None:            
                result = self.__mann_whitney_u(taxon=taxon, group1=group1, group2=group2)
            else:
                result = self.__wilcoxon_signed_rank(taxon=taxon, group1=group1, group2=group2)

            stats_data.append(result)

            significant = result['p-value'] <= self.p_value
            abundant = result['Mean 1 (%)'] >= self.min_abundance_per_group or result['Mean 2 (%)'] >= self.min_abundance_per_group

            if significant and abundant:  # only plot taxa that are significant AND abundant
                
                is_group1_or_group2 = self.taxon_df[GROUP_COLUMN].isin([group1, group2])
                c1 = self.colors[self.groups.index(group1)]
                c2 = self.colors[self.groups.index(group2)]
                p = result['p-value']

                # for paired data, plot all data points, including the ones with missing paired measurements
                Boxplot(self.settings).main(
                    data=self.taxon_df[is_group1_or_group2],
                    x=GROUP_COLUMN,
                    y=taxon,
                    units=self.subject_column,
                    colors=[c1, c2],
                    title=f'{taxon}\n$p = {p:.4f}$',
                    png=f'{dstdir}/{p:.4f}_{taxon}.png'
                )

        self.__save_stats_data(stats_data=stats_data, dstdir=dstdir)

    def __mann_whitney_u(self, taxon: str, group1: str, group2: str) -> Dict[str, Any]:
        is_group1 = self.taxon_df[GROUP_COLUMN] == group1
        is_group2 = self.taxon_df[GROUP_COLUMN] == group2
        x = self.taxon_df.loc[is_group1, taxon]
        y = self.taxon_df.loc[is_group2, taxon]
        statistic, p_value = mannwhitneyu(x=x, y=y)
        return {
            'Taxon': taxon,
            'Group 1': group1,
            'Group 2': group2,
            'Test': 'Mann—Whitney U',
            'N 1': len(x),
            'N 2': len(y),
            'Mean 1 (%)': x.mean(),
            'Mean 2 (%)': y.mean(),
            'p-value': p_value,
        }

    def __wilcoxon_signed_rank(self, taxon: str, group1: str, group2: str) -> Dict[str, Any]:
        # keep only group1 and group2 samples
        in_group1_or_group2 = self.taxon_df[GROUP_COLUMN].isin([group1, group2])
        df = self.taxon_df.loc[in_group1_or_group2, :]

        # pdf: paired data frame, subject x group
        pdf = pd.pivot_table(
            df,
            index=self.subject_column,
            columns=GROUP_COLUMN,
            values=taxon)

        pdf.dropna(inplace=True)  # drop any subjects with missing data, i.e. no paired measurements

        # because there are only group1 and group2 samples, there should be 2 columns for the 2 groups
        # however, if there are less than 2 columns, it means that there are no paired measurements
        n_cols = len(pdf.columns)
        if n_cols < 2:
            n_paired = 0
            p_value = np.nan
        else:
            n_paired = len(pdf)
            x = pdf[group1]
            y = pdf[group2]
            if n_paired == 0:
                p_value = np.nan
            elif (y == 0).all() and (x == 0).all():  # all paired measurements are 0, no change at all
                p_value = np.nan
            else:
                p_value = wilcoxon(x=x, y=y).pvalue

        mean1 = pdf[group1].mean() if group1 in pdf.columns else np.nan
        mean2 = pdf[group2].mean() if group2 in pdf.columns else np.nan
        return {
            'Taxon': taxon,
            'Group 1': group1,
            'Group 2': group2,
            'Test': 'Wilcoxon signed-rank',
            'N Paired': n_paired,
            'Mean 1 (%)': mean1,
            'Mean 2 (%)': mean2,
            'p-value': p_value,
        }

    def __save_stats_data(self, stats_data: List[Dict[str, Any]], dstdir: str):
        stats_df = pd.DataFrame(stats_data).sort_values(
            by='p-value',
            ascending=True
        )

        # only consider abundant taxa for FDR correction
        a = stats_df['Mean 1 (%)'] >= self.min_abundance_per_group
        b = stats_df['Mean 2 (%)'] >= self.min_abundance_per_group
        abundant_taxa = a | b
        p_values = stats_df.loc[abundant_taxa, 'p-value']

        if len(p_values) > 0:
            rejected, pvals_corrected, _, _ = multipletests(
                p_values,
                alpha=0.1,
                method='fdr_bh',  # Benjamini-Hochberg
                is_sorted=False,
                returnsorted=False)
            stats_df.loc[abundant_taxa, 'Benjamini-Hochberg adjusted p-value'] = pvals_corrected
        else:
            stats_df['Benjamini-Hochberg adjusted p-value'] = np.nan  # create an empty column

        stats_df.to_csv(
            f'{dstdir}/statistics.csv',
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
    units: Optional[str]
    colors: List[Tuple[float, float, float, float]]
    title: str
    png: str

    fig: matplotlib.figure.Figure
    ax: matplotlib.axes.Axes

    def main(
            self,
            data: pd.DataFrame,
            x: str,
            y: str,
            units: Optional[str],
            colors: List[Tuple[float, float, float, float]],
            title: str,
            png: str):

        self.data = data
        self.x = x
        self.y = y
        self.units = units
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

        self.fig, self.ax = plt.subplots(figsize=figsize)

    def plot(self):
        sns.boxplot(
            data=self.data,
            x=self.x,
            y=self.y,
            hue=self.x,
            palette=self.colors,
            width=self.BOX_WIDTH,
            linewidth=self.BOX_lINEWIDTH,
            dodge=False,  # to align the boxes on the x axis
            ax=self.ax,
        )
        if self.units is not None:
            sns.lineplot(
                data=self.data,
                x=self.x,
                y=self.y,
                units=self.units,
                estimator=None,
                color='black',
                alpha=0.25,
                linewidth=self.LINEWIDTH,
                sort=False,
                ax=self.ax,
            )
        sns.stripplot(
            data=self.data,
            x=self.x,
            y=self.y,
            hue=self.x,
            palette=self.colors,
            linewidth=self.MARKER_LINEWIDTH,
            ax=self.ax,
        )

    def config(self):
        self.ax.set_title(self.title)
        self.ax.set(xlabel=self.XLABEL, ylabel=self.YLABEL)
        self.ax.xaxis.set_tick_params(width=self.LINEWIDTH)
        self.ax.yaxis.set_tick_params(width=self.LINEWIDTH)
        self.ax.set_ylim(self.YLIM)
        legend = self.ax.get_legend()
        if legend is not None:
            legend.remove()

    def save(self):
        self.fig.tight_layout()
        self.fig.savefig(self.png, dpi=self.DPI)
        plt.close(self.fig)
