import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from sklearn import manifold
from typing import List, Callable, Tuple
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from .tools import edit_fpath
from .template import Processor, Settings


class Ordination(Processor):

    NAME: str
    XY_COLUMNS: Tuple[str, str]

    distance_matrix_tsv: str
    distance_matrix: pd.DataFrame
    sample_coordinate_df: pd.DataFrame
    dstdir: str

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.scatterplot = ScatterPlot(self.settings).main

    def load_distance_matrix(self):
        self.distance_matrix = pd.read_csv(
            self.distance_matrix_tsv,
            sep='\t',
            index_col=0)

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def write_sample_coordinate(self):
        tsv = self.__get_sample_coordinate_fpath(ext='tsv')
        self.sample_coordinate_df.to_csv(tsv, sep='\t')

    def plot_sample_coordinate(self):
        png = self.__get_sample_coordinate_fpath(ext='png')
        self.scatterplot(
            sample_coordinate_df=self.sample_coordinate_df,
            x_column=self.XY_COLUMNS[0],
            y_column=self.XY_COLUMNS[1],
            output_png=png
        )

    def __get_sample_coordinate_fpath(self, ext: str) -> str:
        name = self.NAME.lower().replace('-', '')
        return edit_fpath(
            fpath=self.distance_matrix_tsv,
            old_suffix='.tsv',
            new_suffix=f'-{name}-sample-coordinate.{ext}',
            dstdir=self.dstdir
        )


class PCoA(Ordination):

    NAME = 'PCoA'
    XY_COLUMNS = ('PC1', 'PC2')

    proportion_explained_serise: pd.Series
    proportion_explained_tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, distance_matrix_tsv: str):
        self.distance_matrix_tsv = distance_matrix_tsv

        self.load_distance_matrix()
        self.run_pcoa()
        self.make_dstdir()
        self.write_sample_coordinate()
        self.write_proportion_explained()
        self.plot_sample_coordinate()

    def run_pcoa(self):
        df = self.distance_matrix
        dist_mat = DistanceMatrix(df, list(df.columns))
        result = pcoa(distance_matrix=dist_mat)
        self.sample_coordinate_df = result.samples
        self.proportion_explained_serise = result.proportion_explained

    def write_proportion_explained(self):
        self.proportion_explained_tsv = edit_fpath(
            fpath=self.distance_matrix_tsv,
            old_suffix='.tsv',
            new_suffix=f'-{self.NAME.lower()}-proportion-explained.tsv',
            dstdir=self.dstdir
        )
        self.proportion_explained_serise.to_csv(
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


class NMDS(Ordination):

    NAME = 'NMDS'
    XY_COLUMNS = ('NMDS 1', 'NMDS 2')
    N_COMPONENTS = 2
    METRIC = False  # i.e. non-metric MDS
    N_INIT = 10  # number of independent fitting
    RANDOM_STATE = 1  # to ensure reproducible result
    DISSIMILARITY = 'precomputed'  # distance matrix is precomputed

    embedding: manifold.MDS

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, distance_matrix_tsv: str):
        self.distance_matrix_tsv = distance_matrix_tsv

        self.load_distance_matrix()
        self.run_nmds()
        self.make_dstdir()
        self.write_sample_coordinate()
        self.write_stress()
        self.plot_sample_coordinate()

    def run_nmds(self):
        self.embedding = manifold.MDS(
            n_components=self.N_COMPONENTS,
            metric=self.METRIC,
            n_init=self.N_INIT,
            max_iter=300,
            verbose=0,
            eps=0.001,
            n_jobs=self.threads,
            random_state=self.RANDOM_STATE,
            dissimilarity=self.DISSIMILARITY)

        transformed = self.embedding.fit_transform(
            self.distance_matrix.to_numpy()
        )

        sample_names = list(self.distance_matrix.columns)
        self.sample_coordinate_df = pd.DataFrame(
            data=transformed,
            columns=self.XY_COLUMNS,
            index=sample_names)

    def write_stress(self):
        txt = edit_fpath(
            fpath=self.distance_matrix_tsv,
            old_suffix='.tsv',
            new_suffix='-nmds-stress.txt',
            dstdir=self.dstdir
        )
        with open(txt, 'w') as fh:
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


class TSNE(Ordination):

    NAME = 't-SNE'
    XY_COLUMNS = ('t-SNE 1', 't-SNE 2')
    N_COMPONENTS = 2
    DISTANCE_METRIC = 'precomputed'  # distance matrix is precomputed
    RANDOM_STATE = 1  # to ensure reproducible result
    TSNE_INIT = 'random'  # cannot use PCA becuase distance matrix is precomputed

    embedding: manifold.TSNE

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, distance_matrix_tsv: str):
        self.distance_matrix_tsv = distance_matrix_tsv

        self.load_distance_matrix()
        self.run_tsne()
        self.make_dstdir()
        self.write_sample_coordinate()
        self.plot_sample_coordinate()

    def run_tsne(self):
        self.embedding = manifold.TSNE(
            n_components=self.N_COMPONENTS,
            perplexity=30.0,
            early_exaggeration=12.0,
            learning_rate=200.0,
            n_iter=1000,
            n_iter_without_progress=300,
            min_grad_norm=1e-7,
            metric=self.DISTANCE_METRIC,
            init=self.TSNE_INIT,
            verbose=0,
            random_state=self.RANDOM_STATE,
            method='barnes_hut',
            angle=0.5,
            n_jobs=self.threads,
            square_distances='legacy'
        )

        transformed = self.embedding.fit_transform(
            self.distance_matrix.to_numpy()
        )

        sample_names = list(self.distance_matrix.columns)
        self.sample_coordinate_df = pd.DataFrame(
            data=transformed,
            columns=self.XY_COLUMNS,
            index=sample_names)


class BatchTSNE(Processor):

    distance_matrix_tsvs: List[str]

    tsne: Callable

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.tsne = TSNE(self.settings).main

    def main(self, distance_matrix_tsvs: List[str]):
        self.distance_matrix_tsvs = distance_matrix_tsvs
        for tsv in self.distance_matrix_tsvs:
            self.logger.debug(f't-SNE for {tsv}')
            self.tsne(tsv)


class ScatterPlot(Processor):

    FIGSIZE = (8, 8)
    DPI = 300

    sample_coordinate_df: pd.DataFrame
    x_column: str
    y_column: str
    output_png: str

    ax: Axes

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            sample_coordinate_df: pd.DataFrame,
            x_column: str,
            y_column: str,
            output_png: str):

        self.sample_coordinate_df = sample_coordinate_df
        self.x_column = x_column
        self.y_column = y_column
        self.output_png = output_png

        self.init_figure()
        self.scatterplot()
        self.label_points()
        self.save_figure()

    def init_figure(self):
        plt.figure(figsize=self.FIGSIZE, dpi=self.DPI)

    def scatterplot(self):
        self.ax = sns.scatterplot(
            data=self.sample_coordinate_df,
            x=self.x_column,
            y=self.y_column
        )

    def label_points(self):
        df = self.sample_coordinate_df
        for sample_name in df.index:
            self.ax.text(
                x=df.loc[sample_name, self.x_column],
                y=df.loc[sample_name, self.y_column],
                s=sample_name
            )

    def save_figure(self):
        plt.savefig(self.output_png)
        plt.close()
