from os.path import basename
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

    def main(
            self,
            feature_table_qza: str,
            rooted_tree_qza: str):

        self.feature_table_qza = feature_table_qza
        self.rooted_tree_qza = rooted_tree_qza

        for metric in BETA_METRICS:
            self.run_one_beta_metric_to_tsv(metric=metric)

        for metric in BETA_PHYLOGENETIC_METRICS:
            try:
                self.run_one_beta_phylogenetic_metric_to_tsv(metric=metric)
            except Exception as e:
                self.log_error(metric=metric, exception_instance=e)

    def run_one_beta_metric_to_tsv(self, metric: str):
        RunOneBetaMetric(self.settings).main(
            feature_table_qza=self.feature_table_qza,
            metric=metric)

    def run_one_beta_phylogenetic_metric_to_tsv(self, metric: str):
        RunOneBetaPhylogeneticMetric(self.settings).main(
            feature_table_qza=self.feature_table_qza,
            rooted_tree_qza=self.rooted_tree_qza,
            metric=metric)

    def log_error(self, metric: str, exception_instance: Exception):
        msg = f'"{metric}" results error:\n{exception_instance}'
        self.logger.info(msg)


class RunOneBetaMetric(Processor):

    feature_table_qza: str
    metric: str

    output_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_table_qza: str, metric: str):

        self.feature_table_qza = feature_table_qza
        self.metric = metric

        self.execute()
        self.export()

    def execute(self):
        self.output_qza = f'{self.workdir}/beta-{self.metric}.qza'
        cmd = self.CMD_LINEBREAK.join([
            'qiime diversity beta',
            f'--i-table {self.feature_table_qza}',
            f'--p-metric {self.metric}',
            f'--o-distance-matrix {self.output_qza}'
        ])
        self.call(cmd)

    def export(self):
        ExportBetaDiversityQza(self.settings).main(qza=self.output_qza)


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
            metric: str):

        self.feature_table_qza = feature_table_qza
        self.rooted_tree_qza = rooted_tree_qza
        self.metric = metric

        self.execute()
        self.export()

    def execute(self):
        self.output_qza = f'{self.workdir}/beta-phylogenetic-{self.metric}.qza'
        lines = [
            'qiime diversity beta-phylogenetic',
            f'--i-table {self.feature_table_qza}',
            f'--i-phylogeny {self.rooted_tree_qza}',
            f'--p-metric {self.metric}',
            f'--o-distance-matrix {self.output_qza}'
        ]
        self.call(' \\\n  '.join(lines))

    def export(self):
        ExportBetaDiversityQza(self.settings).main(qza=self.output_qza)


class ExportBetaDiversityQza(Processor):

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
        fname = basename(self.qza)[:-len('.qza')] + '.tsv'
        self.call(f'mv {self.workdir}/distance-matrix.tsv {self.outdir}/{fname}')
