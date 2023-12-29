import os
import pandas as pd
from abc import ABC
from typing import List
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from .tools import edit_fpath
from .template import Processor, Settings
from .exporting import ExportBetaDiversity
from .embedding_core_process import EmbeddingProcessTemplate, TSNECore


class QiimeBetaDiversity(Processor):

    feature_table_qza: str
    rooted_tree_qza: str
    sample_sheet: str
    colormap: str

    distance_matrix_tsvs: List[str]

    def main(
            self,
            feature_table_qza: str,
            rooted_tree_qza: str,
            sample_sheet: str,
            colormap: str):

        self.feature_table_qza = feature_table_qza
        self.rooted_tree_qza = rooted_tree_qza
        self.sample_sheet = sample_sheet
        self.colormap = colormap

        self.run_all_beta_metrics_to_distance_matrix_tsvs()
        self.run_batch_embedding_processes()

    def run_all_beta_metrics_to_distance_matrix_tsvs(self):
        self.distance_matrix_tsvs = RunAllBetaMetricsToDistanceMatrixTsvs(self.settings).main(
            feature_table_qza=self.feature_table_qza,
            rooted_tree_qza=self.rooted_tree_qza)

    def run_batch_embedding_processes(self):
        for Batch in [BatchPCoAProcess, BatchTSNEProcess]:
            Batch(self.settings).main(
                distance_matrix_tsvs=self.distance_matrix_tsvs,
                sample_sheet=self.sample_sheet,
                colormap=self.colormap)


#


class RunAllBetaMetricsToDistanceMatrixTsvs(Processor):

    METRICS = [
        'jaccard',
        'euclidean',
        'braycurtis',
    ]
    PHYLOGENETIC_METRICS = [
        'weighted_unifrac',
        'weighted_normalized_unifrac',
        'generalized_unifrac',
        'unweighted_unifrac'
    ]

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

        for metric in self.METRICS:
            self.run_one_beta_metric_to_tsv(metric=metric)

        for metric in self.PHYLOGENETIC_METRICS:
            self.run_one_beta_phylogenetic_metric_to_tsv(metric=metric)

        self.move_distance_matrix_tsvs_to_outdir()

        return self.distance_matrix_tsvs

    def run_one_beta_metric_to_tsv(self, metric: str):
        try:
            tsv = RunOneBetaMetricToTsv(self.settings).main(
                feature_table_qza=self.feature_table_qza,
                metric=metric)
            self.distance_matrix_tsvs.append(tsv)
        except Exception as e:
            self.log_error(metric=metric, exception_instance=e)

    def run_one_beta_phylogenetic_metric_to_tsv(self, metric: str):
        try:
            tsv = RunOneBetaPhylogeneticMetricToTsv(self.settings).main(
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


class RunOneBetaMetricToTsv(Processor):

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
        log = f'{self.outdir}/qiime-diversity-beta.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime diversity beta',
            f'--i-table {self.feature_table_qza}',
            f'--p-metric {self.metric}',
            f'--o-distance-matrix {self.distance_matrix_qza}',
            f'1>> {log}',
            f'2>> {log}'
        ])
        self.call(cmd)

    def export(self):
        self.distance_matrix_tsv = ExportBetaDiversity(self.settings).main(
            distance_matrix_qza=self.distance_matrix_qza)


class RunOneBetaPhylogeneticMetricToTsv(Processor):

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
        log = f'{self.outdir}/qiime-diversity-beta-phylogenetic.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime diversity beta-phylogenetic',
            f'--i-table {self.feature_table_qza}',
            f'--i-phylogeny {self.rooted_tree_qza}',
            f'--p-metric {self.metric}',
            f'--o-distance-matrix {self.distance_matrix_qza}',
            f'1>> {log}',
            f'2>> {log}'
        ])
        self.call(cmd)

    def export(self):
        self.distance_matrix_tsv = ExportBetaDiversity(self.settings).main(
            distance_matrix_qza=self.distance_matrix_qza)


#


class EmbeddingProcess(EmbeddingProcessTemplate, ABC):

    DSTDIR_NAME = 'beta-embedding'

    def preprocessing(self):
        pass


class PCoAProcess(EmbeddingProcess):

    NAME = 'PCoA'
    XY_COLUMNS = ('PC1', 'PC2')

    proportion_explained_serise: pd.Series

    def main(
            self,
            tsv: str,
            sample_sheet: str,
            colormap: str):

        self.tsv = tsv
        self.sample_sheet = sample_sheet
        self.colormap = colormap

        self.run_main_workflow()
        self.write_proportion_explained()

    def embedding(self):
        df = self.df
        dist_mat = DistanceMatrix(df, list(df.columns))
        result = pcoa(distance_matrix=dist_mat)
        self.sample_coordinate_df = result.samples
        self.proportion_explained_serise = result.proportion_explained

    def write_proportion_explained(self):
        tsv = edit_fpath(
            fpath=self.tsv,
            old_suffix='.tsv',
            new_suffix=f'-{self.NAME.lower()}-proportion-explained.tsv',
            dstdir=self.dstdir
        )
        self.proportion_explained_serise.to_csv(
            tsv,
            sep='\t',
            header=['Proportion Explained']
        )


class TSNEProcess(EmbeddingProcess):

    NAME = 't-SNE'
    XY_COLUMNS = ('t-SNE 1', 't-SNE 2')

    def main(
            self,
            tsv: str,
            sample_sheet: str,
            colormap: str):

        self.tsv = tsv
        self.sample_sheet = sample_sheet
        self.colormap = colormap

        self.run_main_workflow()

    def embedding(self):
        self.sample_coordinate_df = TSNECore(self.settings).main(
            df=self.df,
            data_structure='distance_matrix'
        )


#


class BatchEmbeddingProcess(Processor):

    distance_matrix_tsvs: List[str]
    sample_sheet: str
    colormap: str

    embedding: EmbeddingProcessTemplate

    def main(
            self,
            distance_matrix_tsvs: List[str],
            sample_sheet: str,
            colormap: str):

        self.distance_matrix_tsvs = distance_matrix_tsvs
        self.sample_sheet = sample_sheet
        self.colormap = colormap

        for tsv in self.distance_matrix_tsvs:
            self.logger.debug(f'{self.embedding.NAME} for {tsv}')
            self.embedding.main(
                tsv=tsv,
                sample_sheet=self.sample_sheet,
                colormap=self.colormap)


class BatchPCoAProcess(BatchEmbeddingProcess):

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.embedding = PCoAProcess(self.settings)


class BatchTSNEProcess(BatchEmbeddingProcess):

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.embedding = TSNEProcess(self.settings)
