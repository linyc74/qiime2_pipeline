import os
from typing import List
from .tools import edit_fpath
from .template import Processor
from .exporting import ExportBetaDiversity


BETA_METRICS = [
    'jaccard',
    'euclidean',
    'braycurtis',
    'cosine',
    'correlation',
]
BETA_PHYLOGENETIC_METRICS = [
    'weighted_unifrac',
    'weighted_normalized_unifrac',
    'generalized_unifrac',
    'unweighted_unifrac'
]


class BetaDiversity(Processor):

    feature_table_qza: str
    rooted_tree_qza: str

    distance_matrix_tsvs: List[str]

    def main(
            self,
            feature_table_qza: str,
            rooted_tree_qza: str):

        self.feature_table_qza = feature_table_qza
        self.rooted_tree_qza = rooted_tree_qza

        self.distance_matrix_tsvs = []

        for metric in BETA_METRICS:
            self.run_one_beta_metric_to_tsv(metric=metric)

        for metric in BETA_PHYLOGENETIC_METRICS:
            self.run_one_beta_phylogenetic_metric_to_tsv(metric=metric)

        self.move_distance_matrix_tsvs_to_outdir()

        return self.distance_matrix_tsvs

    def run_one_beta_metric_to_tsv(self, metric: str):
        tsv = RunOneBetaMetric(self.settings).main(
            feature_table_qza=self.feature_table_qza,
            metric=metric)
        self.distance_matrix_tsvs.append(tsv)

    def run_one_beta_phylogenetic_metric_to_tsv(self, metric: str):
        try:
            tsv = RunOneBetaPhylogeneticMetric(self.settings).main(
                feature_table_qza=self.feature_table_qza,
                rooted_tree_qza=self.rooted_tree_qza,
                metric=metric)
            self.distance_matrix_tsvs.append(tsv)
        except Exception as e:
            self.log_error(metric=metric, exception_instance=e)

    def move_distance_matrix_tsvs_to_outdir(self):
        dstdir = f'{self.outdir}/beta-diversity'
        os.makedirs(dstdir, exist_ok=True)
        for i, tsv in enumerate(self.distance_matrix_tsvs):
            new = edit_fpath(
                fpath=tsv,
                old_suffix='',
                new_suffix='',
                dstdir=dstdir)
            self.call(f'mv {tsv} {new}')
            self.distance_matrix_tsvs[i] = new

    def log_error(self, metric: str, exception_instance: Exception):
        msg = f'"{metric}" results error:\n{exception_instance}'
        self.logger.info(msg)


class RunOneBetaMetric(Processor):

    feature_table_qza: str
    metric: str

    distance_matrix_qza: str
    distance_matrix_tsv: str

    def main(self, feature_table_qza: str, metric: str) -> str:
        self.feature_table_qza = feature_table_qza
        self.metric = metric
        self.execute()
        self.export()
        return self.distance_matrix_tsv

    def execute(self):
        self.distance_matrix_qza = f'{self.workdir}/{self.metric}.qza'
        cmd = self.CMD_LINEBREAK.join([
            'qiime diversity beta',
            f'--i-table {self.feature_table_qza}',
            f'--p-metric {self.metric}',
            f'--o-distance-matrix {self.distance_matrix_qza}'
        ])
        self.call(cmd)

    def export(self):
        self.distance_matrix_tsv = ExportBetaDiversity(self.settings).main(
            distance_matrix_qza=self.distance_matrix_qza)


class RunOneBetaPhylogeneticMetric(Processor):

    feature_table_qza: str
    rooted_tree_qza: str
    metric: str

    distance_matrix_qza: str
    distance_matrix_tsv: str

    def main(
            self,
            feature_table_qza: str,
            rooted_tree_qza: str,
            metric: str) -> str:
        self.feature_table_qza = feature_table_qza
        self.rooted_tree_qza = rooted_tree_qza
        self.metric = metric
        self.execute()
        self.export()
        return self.distance_matrix_tsv

    def execute(self):
        self.distance_matrix_qza = f'{self.workdir}/{self.metric}.qza'
        cmd = self.CMD_LINEBREAK.join([
            'qiime diversity beta-phylogenetic',
            f'--i-table {self.feature_table_qza}',
            f'--i-phylogeny {self.rooted_tree_qza}',
            f'--p-metric {self.metric}',
            f'--o-distance-matrix {self.distance_matrix_qza}'
        ])
        self.call(cmd)

    def export(self):
        self.distance_matrix_tsv = ExportBetaDiversity(self.settings).main(
            distance_matrix_qza=self.distance_matrix_qza)
