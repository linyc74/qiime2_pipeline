import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from sklearn import manifold
from typing import List, Tuple
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from .tools import edit_fpath
from .grouping import AddGroupColumn
from .template import Processor, Settings


class Ordination(Processor):

    NAME: str
    XY_COLUMNS: Tuple[str, str]
    GROUP_COLUMN: str = AddGroupColumn.GROUP_COLUMN

    distance_matrix_tsv: str
    group_keywords: List[str]

    distance_matrix: pd.DataFrame
    sample_coordinate_df: pd.DataFrame
    dstdir: str

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.scatterplot = ScatterPlot(self.settings).main

    def main(
            self,
            distance_matrix_tsv: str,
            group_keywords: List[str]):
        pass

    def run_main_workflow(self):
        self.load_distance_matrix()
        self.run_dim_reduction()
        self.add_group_column()
        self.make_dstdir()
        self.write_sample_coordinate()
        self.plot_sample_coordinate()

    def load_distance_matrix(self):
        self.distance_matrix = pd.read_csv(
            self.distance_matrix_tsv,
            sep='\t',
            index_col=0)

    def run_dim_reduction(self):
        pass

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def add_group_column(self):
        self.sample_coordinate_df = AddGroupColumn(self.settings).main(
            df=self.sample_coordinate_df,
            group_keywords=self.group_keywords)

    def write_sample_coordinate(self):
        tsv = self.__get_sample_coordinate_fpath(ext='tsv')
        self.sample_coordinate_df.to_csv(tsv, sep='\t')

    def plot_sample_coordinate(self):
        png = self.__get_sample_coordinate_fpath(ext='png')
        self.scatterplot(
            sample_coordinate_df=self.sample_coordinate_df,
            x_column=self.XY_COLUMNS[0],
            y_column=self.XY_COLUMNS[1],
            hue_column=self.GROUP_COLUMN,
            output_png=png)

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

    def main(
            self,
            distance_matrix_tsv: str,
            group_keywords: List[str]):

        self.distance_matrix_tsv = distance_matrix_tsv
        self.group_keywords = group_keywords

        self.run_main_workflow()
        self.write_proportion_explained()

    def run_dim_reduction(self):
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


class NMDS(Ordination):

    NAME = 'NMDS'
    XY_COLUMNS = ('NMDS 1', 'NMDS 2')
    N_COMPONENTS = 2
    METRIC = False  # i.e. non-metric MDS
    N_INIT = 10  # number of independent fitting
    EPS = 1e-3  # stress tolerance for convergence
    RANDOM_STATE = 1  # to ensure reproducible result
    DISSIMILARITY = 'precomputed'  # distance matrix is precomputed

    embedding: manifold.MDS
    stress: float

    def main(
            self,
            distance_matrix_tsv: str,
            group_keywords: List[str]):

        self.distance_matrix_tsv = distance_matrix_tsv
        self.group_keywords = group_keywords

        self.run_main_workflow()
        self.normalize_stress()
        self.write_stress()

    def run_dim_reduction(self):
        self.embedding = manifold.MDS(
            n_components=self.N_COMPONENTS,
            metric=self.METRIC,
            n_init=self.N_INIT,
            max_iter=300,
            verbose=0,
            eps=self.EPS,
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

    def normalize_stress(self):
        """
        Normalize with sum of squared distances
        i.e. Kruskal Stress, or Stress_1
        https://stackoverflow.com/questions/36428205/stress-attribute-sklearn-manifold-mds-python
        """
        distances = self.distance_matrix.to_numpy()
        squared_sum = np.sum(distances ** 2) / 2  # diagonal symmetry, thus divide by 2
        self.stress = np.sqrt(self.embedding.stress_ / squared_sum)

    def write_stress(self):
        txt = edit_fpath(
            fpath=self.distance_matrix_tsv,
            old_suffix='.tsv',
            new_suffix='-nmds-stress.txt',
            dstdir=self.dstdir)
        with open(txt, 'w') as fh:
            fh.write(str(self.stress))


class TSNE(Ordination):

    NAME = 't-SNE'
    XY_COLUMNS = ('t-SNE 1', 't-SNE 2')
    N_COMPONENTS = 2
    PERPLEXITY = 3.0
    DISTANCE_METRIC = 'precomputed'  # distance matrix is precomputed
    RANDOM_STATE = 1  # to ensure reproducible result
    TSNE_INIT = 'random'  # cannot use PCA becuase distance matrix is precomputed

    embedding: manifold.TSNE

    def main(
            self,
            distance_matrix_tsv: str,
            group_keywords: List[str]):

        self.distance_matrix_tsv = distance_matrix_tsv
        self.group_keywords = group_keywords

        self.run_main_workflow()

    def run_dim_reduction(self):
        self.embedding = manifold.TSNE(
            n_components=self.N_COMPONENTS,
            perplexity=self.PERPLEXITY,
            early_exaggeration=12.0,
            learning_rate=200.0,
            n_iter=1000,
            n_iter_without_progress=300,
            min_grad_norm=1e-7,
            metric=self.DISTANCE_METRIC,
            init=self.TSNE_INIT,
            verbose=1,
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


class BatchOrdination(Processor):

    distance_matrix_tsvs: List[str]
    group_keywords: List[str]

    ordination: Ordination

    def main(
            self,
            distance_matrix_tsvs: List[str],
            group_keywords: List[str]):

        self.distance_matrix_tsvs = distance_matrix_tsvs
        self.group_keywords = group_keywords

        for tsv in self.distance_matrix_tsvs:
            self.logger.debug(f'{self.ordination.NAME} for {tsv}')
            self.ordination.main(
                distance_matrix_tsv=tsv,
                group_keywords=self.group_keywords)


class BatchPCoA(BatchOrdination):

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.ordination = PCoA(self.settings)


class BatchNMDS(BatchOrdination):

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.ordination = NMDS(self.settings)


class BatchTSNE(BatchOrdination):

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.ordination = TSNE(self.settings)


class ScatterPlot(Processor):

    FIGSIZE = (8, 8)
    DPI = 300

    sample_coordinate_df: pd.DataFrame
    x_column: str
    y_column: str
    group_column: str
    output_png: str

    ax: Axes

    def main(
            self,
            sample_coordinate_df: pd.DataFrame,
            x_column: str,
            y_column: str,
            hue_column: str,
            output_png: str):

        self.sample_coordinate_df = sample_coordinate_df
        self.x_column = x_column
        self.y_column = y_column
        self.group_column = hue_column
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
            y=self.y_column,
            hue=self.group_column)

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
