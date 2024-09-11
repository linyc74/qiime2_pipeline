import pandas as pd
from .template import Processor
from .exporting import ExportFeatureTable


class AlphaRarefaction(Processor):

    METRICS = ['observed_features', 'shannon', 'simpson']
    MIN_DEPTH = 1
    STEPS = 20
    ITERATIONS = 10

    feature_table_qza: str

    max_depth: int

    def main(self, feature_table_qza: str):
        self.feature_table_qza = feature_table_qza
        self.set_max_depth()
        self.run_rarefaction()

    def set_max_depth(self):
        tsv = ExportFeatureTable(self.settings).main(feature_table_qza=self.feature_table_qza)
        df = pd.read_csv(tsv, sep='\t', index_col=0)
        sum_per_column = df.sum()
        self.max_depth = int(max(sum_per_column) * 0.99)  # 99%, slightly lower than the max depth

    def run_rarefaction(self):
        metrics_line = ' '.join([f'--p-metrics {m}' for m in self.METRICS])
        output_qzv = f'{self.outdir}/rarefaction.qzv'
        log = f'{self.outdir}/qiime-diversity-alpha-rarefaction.log'
        args = [
            'qiime diversity alpha-rarefaction',
            f'--i-table "{self.feature_table_qza}"',
            metrics_line,
            f'--p-max-depth {self.max_depth}',
            f'--p-min-depth {self.MIN_DEPTH}',
            f'--p-steps {self.STEPS}',
            f'--p-iterations {self.ITERATIONS}',
            f'--o-visualization "{output_qzv}"',
            f'1> "{log}"',
            f'2> "{log}"'
        ]
        self.call(self.CMD_LINEBREAK.join(args))
