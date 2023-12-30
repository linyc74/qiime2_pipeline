import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List
from .template import Processor
from .grouping import AddGroupColumn, GROUP_COLUMN


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
    ALPHA_DIVERSITY_DIRNAME = 'alpha-diversity'

    feature_table_qza: str
    sample_sheet: str
    alpha_metrics: List[str]
    colors: list

    df: pd.DataFrame

    def main(
            self,
            feature_table_qza: str,
            sample_sheet: str,
            alpha_metrics: List[str],
            colors: list):

        self.feature_table_qza = feature_table_qza
        self.sample_sheet = sample_sheet
        self.alpha_metrics = self.ALPHA_METRICS if alpha_metrics == [] else alpha_metrics
        self.colors = colors

        self.df = pd.DataFrame()
        for metric in self.alpha_metrics:
            self.run_one(metric=metric)
        self.add_group_column()
        self.save_csv()
        self.plot()

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

    def add_group_column(self):
        self.df = AddGroupColumn(self.settings).main(
            df=self.df,
            sample_sheet=self.sample_sheet)

    def save_csv(self):
        dstdir = f'{self.outdir}/{self.ALPHA_DIVERSITY_DIRNAME}'
        os.makedirs(dstdir, exist_ok=True)
        self.df.to_csv(f'{dstdir}/alpha-diversity.csv', index=True)

    def plot(self):
        PlotAlphaDiversity(self.settings).main(
            df=self.df,
            dstdir=f'{self.outdir}/{self.ALPHA_DIVERSITY_DIRNAME}',
            colors=self.colors)


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
            f'1>> {log}',
            f'2>> {log}'
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
            f'1>> {log}',
            f'2>> {log}'
        ])
        self.call(cmd)


class PlotAlphaDiversity(Processor):

    FIGSIZE = (10 / 2.54, 10 / 2.54)
    BOX_WIDTH = 0.5
    DPI = 600

    df: pd.DataFrame
    dstdir: str
    colors: list

    alpha_metrics: List[str]

    def main(
            self,
            df: pd.DataFrame,
            dstdir: str,
            colors: list):

        self.df = df
        self.dstdir = dstdir
        self.colors = colors

        self.set_alpha_metrics()
        for metric in self.alpha_metrics:
            self.plot_one(metric=metric)

    def set_alpha_metrics(self):
        self.alpha_metrics = [
            c for c in self.df.columns if c != GROUP_COLUMN
        ]

    def plot_one(self, metric: str):

        plt.figure(figsize=self.FIGSIZE, dpi=self.DPI)
        sns.boxplot(
            data=self.df,
            x=GROUP_COLUMN,
            y=metric,
            palette=self.colors,
            width=self.BOX_WIDTH
        )
        for ext in ['pdf', 'png']:
            plt.savefig(f'{self.dstdir}/{metric}.{ext}', dpi=self.DPI)
        plt.close()
