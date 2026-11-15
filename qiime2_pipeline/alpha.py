import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import mannwhitneyu, wilcoxon
from typing import List, Tuple, Dict, Union, Optional
from .template import Processor
from .grouping import AddGroupColumn, GROUP_COLUMN


DSTDIR_NAME = 'alpha-diversity'


class AlphaDiversity(Processor):

    ALPHA_METRICS = [
        'chao1',
        'shannon',
        'gini_index',
        'mcintosh_e',
        'pielou_e',
        'simpson',
        'observed_features',
        'fisher_alpha',
    ]

    feature_table_qza: str
    sample_sheet: str
    subject_column: Optional[str]
    alpha_metrics: List[str]
    colors: list

    df: pd.DataFrame

    def main(
            self,
            feature_table_qza: str,
            sample_sheet: str,
            subject_column: Optional[str],
            alpha_metrics: List[str],
            colors: list):

        self.feature_table_qza = feature_table_qza
        self.sample_sheet = sample_sheet
        self.subject_column = subject_column
        self.alpha_metrics = self.ALPHA_METRICS if alpha_metrics == [] else alpha_metrics
        self.colors = colors

        self.df = pd.DataFrame()
        for metric in self.alpha_metrics:
            self.run_one(metric=metric)

        self.df = AddGroupColumn(self.settings).main(df=self.df, sample_sheet=self.sample_sheet)

        if self.subject_column is not None:
            self.add_subject_column()

        self.reorder_samples_by_sample_sheet()

        os.makedirs(f'{self.outdir}/{DSTDIR_NAME}', exist_ok=True)
        self.df.to_csv(f'{self.outdir}/{DSTDIR_NAME}/alpha-diversity.csv', index=True)

        Plot(self.settings).main(df=self.df, colors=self.colors, subject_column=self.subject_column)

        Statistics(self.settings).main(df=self.df, subject_column=self.subject_column)

    def run_one(self, metric: str):
        qza = RunOneAlphaMetric(self.settings).main(
            feature_table_qza=self.feature_table_qza,
            metric=metric)

        df = ReadAlphaDiversityQza(self.settings).main(qza=qza)

        self.df = self.df.merge(
            right=df,
            how='outer',
            left_index=True,
            right_index=True)

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

    def reorder_samples_by_sample_sheet(self):
        sample_order = pd.read_csv(self.sample_sheet, index_col=0).index
        self.df = self.df.loc[sample_order, :]


class RunOneAlphaMetric(Processor):

    feature_table_qza: str
    metric: str
    output_qza: str

    def main(
            self,
            feature_table_qza: str,
            metric: str) -> str:

        self.feature_table_qza = feature_table_qza
        self.metric = metric

        self.set_output_qza()
        self.execute()

        return self.output_qza

    def set_output_qza(self):
        self.output_qza = f'{self.workdir}/alpha-{self.metric}.qza'

    def execute(self):
        log = f'{self.outdir}/qiime-diversity-alpha.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime diversity alpha',
            f'--i-table {self.feature_table_qza}',
            f'--p-metric {self.metric}',
            f'--o-alpha-diversity {self.output_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)


class ReadAlphaDiversityQza(Processor):

    qza: str

    def main(self, qza: str) -> pd.DataFrame:

        self.qza = qza
        self.qza_to_tsv()
        return pd.read_csv(
            f'{self.workdir}/alpha-diversity.tsv',
            sep='\t',
            index_col=0
        )

    def qza_to_tsv(self):
        log = f'{self.outdir}/qiime-tools-export.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools export',
            f'--input-path {self.qza}',
            f'--output-path {self.workdir}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)


class Plot(Processor):

    df: pd.DataFrame
    colors: list
    subject_column: Optional[str]

    figsize: Tuple[float, float]
    box_width: float
    dpi: int
    line_width: float
    fontsize: int

    alpha_metrics: List[str]

    def main(self, df: pd.DataFrame, colors: list, subject_column: Optional[str]):
        self.df = df
        self.colors = colors
        self.subject_column = subject_column

        self.set_parameters()
        self.set_alpha_metrics()
        for metric in self.alpha_metrics:
            self.plot_one(metric=metric)

    def set_parameters(self):
        n_groups = self.df[GROUP_COLUMN].nunique()
        if self.settings.for_publication:
            self.figsize = get_figsize_for_publication(n_groups=n_groups)
            self.box_width = 0.35
            self.dpi = 600
            self.line_width = 0.5
            self.fontsize = 7
        else:
            self.figsize = get_figsize(n_groups=n_groups)
            self.box_width = 0.5
            self.dpi = 300
            self.line_width = 1.0
            self.fontsize = 12

    def set_alpha_metrics(self):
        self.alpha_metrics = [
            c for c in self.df.columns if c not in [GROUP_COLUMN, self.subject_column]
        ]

    def plot_one(self, metric: str):
        plt.rcParams['font.size'] = self.fontsize
        plt.rcParams['axes.linewidth'] = self.line_width

        plt.figure(figsize=self.figsize, dpi=self.dpi)

        sns.boxplot(
            data=self.df,
            x=GROUP_COLUMN,
            y=metric,
            palette=self.colors,
            width=self.box_width,
            linewidth=self.line_width
        )
        if self.subject_column is not None:
            sns.lineplot(
                data=self.df,
                x=GROUP_COLUMN,
                y=metric,
                units=self.subject_column,
                estimator=None,
                color='black',
                alpha=0.25,
                linewidth=self.line_width,
                sort=False
            )
        sns.stripplot(
            data=self.df,
            x=GROUP_COLUMN,
            y=metric,
            palette=self.colors,
            linewidth=self.line_width,
        )

        plt.ylabel(metric.replace('_', ' ').title())  # e.g. 'observed_features' -> 'Observed Features'

        plt.tight_layout()

        os.makedirs(f'{self.outdir}/{DSTDIR_NAME}', exist_ok=True)
        for ext in ['pdf', 'png']:
            plt.savefig(f'{self.outdir}/{DSTDIR_NAME}/{metric}.{ext}', dpi=self.dpi)

        plt.close()


def get_figsize(n_groups: int) -> Tuple[float, float]:
    w = (3 * n_groups + 3) / 2.54
    h = 10 / 2.54
    return w, h


def get_figsize_for_publication(n_groups: int) -> Tuple[float, float]:
    w = (1.2 * n_groups + 2) / 2.54
    h = 4 / 2.54
    return w, h


class Statistics(Processor):

    df: pd.DataFrame
    subject_column: Optional[str]
    stats_data: List[Dict[str, Union[str, float]]]

    def main(self, df: pd.DataFrame, subject_column: Optional[str]):
        self.df = df
        self.subject_column = subject_column

        alpha_metrics = [c for c in self.df.columns if c not in [GROUP_COLUMN, self.subject_column]]
        groups = self.df[GROUP_COLUMN].unique()

        self.stats_data = []

        if self.subject_column is None:
            test = self.mann_whitney_u
        else:
            test = self.wilcoxon_signed_rank

        for metric in alpha_metrics:
            for group1, group2 in combinations(groups, 2):
                test(metric=metric, group1=group1, group2=group2)

        os.makedirs(f'{self.outdir}/{DSTDIR_NAME}', exist_ok=True)     
        pd.DataFrame(self.stats_data).to_csv(
            f'{self.outdir}/{DSTDIR_NAME}/alpha-diversity-statistics.csv',
            index=False)

    def mann_whitney_u(self, metric: str, group1: str, group2: str):
        x = self.df.loc[self.df[GROUP_COLUMN] == group1, metric]
        y = self.df.loc[self.df[GROUP_COLUMN] == group2, metric]
        result = mannwhitneyu(x=x, y=y)
        self.stats_data.append({
            'Diversity': metric,
            'Group 1': group1,
            'Group 2': group2,
            'Test': 'Mann—Whitney U',
            'N 1': len(x),
            'N 2': len(y),
            'Mean 1': x.mean(),
            'Mean 2': y.mean(),
            'p-value': result.pvalue,
        })
    
    def wilcoxon_signed_rank(self, metric: str, group1: str, group2: str):
        # keep only group1 and group2 samples
        in_group1_or_group2 = self.df[GROUP_COLUMN].isin([group1, group2])
        df = self.df.loc[in_group1_or_group2, :]

        # pdf: paired data frame, subject x group
        pdf = pd.pivot_table(
            df,
            index=self.subject_column,
            columns=GROUP_COLUMN,
            values=metric)

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
        self.stats_data.append({
            'Diversity': metric,
            'Group 1': group1,
            'Group 2': group2,
            'Test': 'Wilcoxon signed-rank',
            'N Paired': n_paired,
            'Mean 1': mean1,
            'Mean 2': mean2,
            'p-value': p_value
        })
