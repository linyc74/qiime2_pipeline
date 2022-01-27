import os
import numpy as np
import pandas as pd
from sklearn import manifold
from typing import List, Callable
from skbio.stats.ordination import pcoa
from skbio import DistanceMatrix, OrdinationResults
from .tools import edit_fpath
from .template import Processor, Settings


class PCoA(Processor):

    distance_matrix_tsv: str

    distance_matrix: DistanceMatrix
    result: OrdinationResults
    pcoa_outdir: str
    sample_coordinate_tsv: str
    proportion_explained_tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, distance_matrix_tsv: str):
        self.distance_matrix_tsv = distance_matrix_tsv

        self.load_distance_matrix()
        self.run_pcoa()
        self.make_pcoa_outdir()
        self.write_sample_coordinate()
        self.write_proportion_explained()

    def load_distance_matrix(self):
        df = pd.read_csv(
            self.distance_matrix_tsv,
            sep='\t',
            index_col=0
        )
        self.distance_matrix = DistanceMatrix(
            df, list(df.columns)
        )

    def run_pcoa(self):
        self.result = pcoa(distance_matrix=self.distance_matrix)

    def make_pcoa_outdir(self):
        self.pcoa_outdir = f'{self.outdir}/PCoA'
        os.makedirs(self.pcoa_outdir, exist_ok=True)

    def write_sample_coordinate(self):
        self.sample_coordinate_tsv = edit_fpath(
            fpath=self.distance_matrix_tsv,
            old_suffix='.tsv',
            new_suffix='-pcoa-sample-coordinate.tsv',
            dstdir=self.pcoa_outdir
        )
        self.result.samples.to_csv(self.sample_coordinate_tsv, sep='\t')

    def write_proportion_explained(self):
        self.proportion_explained_tsv = edit_fpath(
            fpath=self.distance_matrix_tsv,
            old_suffix='.tsv',
            new_suffix='-pcoa-proportion-explained.tsv',
            dstdir=self.pcoa_outdir
        )
        self.result.proportion_explained.to_csv(
            self.proportion_explained_tsv,
            sep='\t',
            header=['Proportion Explained']
        )


class BatchPCoA(Processor):

    distance_matrix_tsvs: List[str]

    pcoa: Callable

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.pcoa = PCoA(self.settings).main

    def main(self, distance_matrix_tsvs: List[str]):
        self.distance_matrix_tsvs = distance_matrix_tsvs
        for tsv in self.distance_matrix_tsvs:
            self.logger.debug(f'PCoA for {tsv}')
            self.pcoa(tsv)


class NMDS(Processor):

    N_COMPONENTS = 2
    METRIC = False  # i.e. non-metric MDS
    N_INIT = 10  # number of independent fitting
    DISSIMILARITY = 'precomputed'  # distance matrix is precomputed
    NMDS_COLUMNS = ['NMDS1', 'NMDS2']

    distance_matrix_tsv: str

    distance_matrix: pd.DataFrame
    embedding: manifold.MDS
    sample_coordinate_df: pd.DataFrame
    nmds_outdir: str
    sample_coordinate_tsv: str
    stress_txt: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, distance_matrix_tsv: str):
        self.distance_matrix_tsv = distance_matrix_tsv

        self.load_distance_matrix()
        self.run_nmds()
        self.make_nmds_outdir()
        self.write_sample_coordinate()
        self.write_stress()

    def load_distance_matrix(self):
        self.distance_matrix = pd.read_csv(
            self.distance_matrix_tsv,
            sep='\t',
            index_col=0
        )

    def run_nmds(self):
        self.embedding = manifold.MDS(
            n_components=self.N_COMPONENTS,
            metric=self.METRIC,
            n_init=self.N_INIT,
            max_iter=300,
            verbose=0,
            eps=0.001,
            n_jobs=self.threads,
            random_state=None,
            dissimilarity=self.DISSIMILARITY)

        transformed = self.embedding.fit_transform(
            self.distance_matrix.to_numpy()
        )

        sample_names = list(self.distance_matrix.columns)
        self.sample_coordinate_df = pd.DataFrame(
            data=transformed,
            columns=self.NMDS_COLUMNS,
            index=sample_names)

    def make_nmds_outdir(self):
        self.nmds_outdir = f'{self.outdir}/NMDS'
        os.makedirs(self.nmds_outdir, exist_ok=True)

    def write_sample_coordinate(self):
        self.sample_coordinate_tsv = edit_fpath(
            fpath=self.distance_matrix_tsv,
            old_suffix='.tsv',
            new_suffix='-nmds-sample-coordinate.tsv',
            dstdir=self.nmds_outdir
        )
        self.sample_coordinate_df.to_csv(self.sample_coordinate_tsv, sep='\t')

    def write_stress(self):
        self.stress_txt = edit_fpath(
            fpath=self.distance_matrix_tsv,
            old_suffix='.tsv',
            new_suffix='-nmds-stress.txt',
            dstdir=self.nmds_outdir
        )
        with open(self.stress_txt, 'w') as fh:
            fh.write(str(self.embedding.stress_))


class BatchNMDS(Processor):

    distance_matrix_tsvs: List[str]

    nmds: Callable

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.nmds = NMDS(self.settings).main

    def main(self, distance_matrix_tsvs: List[str]):
        self.distance_matrix_tsvs = distance_matrix_tsvs
        for tsv in self.distance_matrix_tsvs:
            self.logger.debug(f'NMDS for {tsv}')
            self.nmds(tsv)
