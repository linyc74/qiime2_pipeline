import os
from .template import Processor, Settings


BETA_METRICS = [
    'jaccard', 'euclidean', 'cosine',
]
BETA_PHYLOGENETIC_METRICS = [
    'weighted_normalized_unifrac', 'weighted_unifrac',
    'generalized_unifrac', 'unweighted_unifrac'
]


class BetaDiversity(Processor):

    feature_table_qza: str
    rooted_tree_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.run_one_beta_metric = RunOneBetaMetric(self.settings).main
        self.run_one_beta_phylogenetic_metric = RunOneBetaPhylogeneticMetric(self.settings).main
        self.qza_to_tsv = QzaToTsv(self.settings).main

    def main(
            self,
            feature_tabe_qza: str,
            rooted_tree_qza: str):

        self.feature_table_qza = feature_tabe_qza
        self.rooted_tree_qza = rooted_tree_qza

        for metric in BETA_METRICS:
            qza = self.run_one_beta_metric(
                feature_table_qza=self.feature_table_qza,
                metric=metric)
            self.qza_to_tsv(qza=qza)

        for metric in BETA_PHYLOGENETIC_METRICS:
            qza = self.run_one_beta_phylogenetic_metric(
                feature_table_qza=self.feature_table_qza,
                rooted_tree_qza=self.rooted_tree_qza,
                metric=metric)
            self.qza_to_tsv(qza=qza)


class RunOneBetaMetric(Processor):

    feature_table_qza: str
    metric: str

    output_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_table_qza: str, metric: str) -> str:

        self.feature_table_qza = feature_table_qza
        self.metric = metric

        self.set_output_qza()
        self.execute()

        return self.output_qza

    def set_output_qza(self):
        self.output_qza = f'{self.workdir}/beta-{self.metric}.qza'

    def execute(self):
        lines = [
            'qiime diversity beta',
            f'--i-table {self.feature_table_qza}',
            f'--p-metric {self.metric}',
            f'--o-distance-matrix {self.output_qza}'
        ]
        self.call(' \\\n  '.join(lines))


class RunOneBetaPhylogeneticMetric(Processor):

    feature_table_qza: str
    rooted_tree_qza: str
    metric: str

    output_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            feature_table_qza: str,
            rooted_tree_qza: str,
            metric: str) -> str:

        self.feature_table_qza = feature_table_qza
        self.rooted_tree_qza = rooted_tree_qza
        self.metric = metric

        self.set_output_qza()
        self.execute()

        return self.output_qza

    def set_output_qza(self):
        self.output_qza = f'{self.workdir}/beta-phylogenetic-{self.metric}.qza'

    def execute(self):
        lines = [
            'qiime diversity beta-phylogenetic',
            f'--i-table {self.feature_table_qza}',
            f'--i-phylogeny {self.rooted_tree_qza}',
            f'--p-metric {self.metric}',
            f'--o-distance-matrix {self.output_qza}'
        ]
        self.call(' \\\n  '.join(lines))


class QzaToTsv(Processor):

    qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, qza: str):

        self.qza = qza
        self.qza_to_tsv()
        self.move_tsv()

    def qza_to_tsv(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def move_tsv(self):
        fname = os.path.basename(self.qza)[:-len('.qza')] + '.tsv'
        os.rename(
            f'{self.workdir}/distance-matrix.tsv',
            f'{self.outdir}/{fname}'
        )
