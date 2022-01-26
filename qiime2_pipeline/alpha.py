import pandas as pd
from .template import Processor, Settings


ALPHA_METRICS = [
    'chao1', 'osd', 'singles', 'shannon', 'gini_index', 'mcintosh_e',
    'michaelis_menten_fit', 'chao1_ci', 'pielou_e', 'simpson',
    'observed_features', 'fisher_alpha',
]


class AlphaDiversity(Processor):

    feature_table_qza: str
    df: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            feature_table_qza: str):

        self.feature_table_qza = feature_table_qza
        self.df = pd.DataFrame()

        for metric in ALPHA_METRICS:
            self.run_one(metric)

        self.df.to_csv(f'{self.outdir}/alpha-diversity.csv', index=True)

    def run_one(self, metric: str):
        qza = RunOneAlphaMetric(self.settings).main(
            feature_table_qza=self.feature_table_qza, metric=metric)

        df = ReadAlphaDiversityQza(self.settings).main(qza=qza)

        self.df = self.df.merge(
            right=df,
            how='outer',
            left_index=True,
            right_index=True)


class RunOneAlphaMetric(Processor):

    feature_table_qza: str
    metric: str
    output_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

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
        cmd = self.CMD_LINEBREAK.join([
            'qiime diversity alpha',
            f'--i-table {self.feature_table_qza}',
            f'--p-metric {self.metric}',
            f'--o-alpha-diversity {self.output_qza}'
        ])
        self.call(cmd)


class ReadAlphaDiversityQza(Processor):

    qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, qza: str) -> pd.DataFrame:

        self.qza = qza
        self.qza_to_tsv()
        return pd.read_csv(
            f'{self.workdir}/alpha-diversity.tsv',
            sep='\t',
            index_col=0
        )

    def qza_to_tsv(self):
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools export',
            f'--input-path {self.qza}',
            f'--output-path {self.workdir}',
        ])
        self.call(cmd)
