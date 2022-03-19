import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import manifold
from matplotlib.axes import Axes
from sklearn import decomposition
from abc import ABC, abstractmethod
from typing import List, Tuple, Union
from .template import Processor


class Core(Processor, ABC):

    XY_COLUMNS: Tuple[float, float]
    DATA_STRUCTURES = [
        'distance_matrix',
        'row_features',
        'column_features',
    ]
    N_COMPONENTS = 2
    RANDOM_STATE = 1  # to ensure reproducible result

    df: pd.DataFrame
    data_structure: str

    embedding: Union[manifold.TSNE, manifold.MDS, decomposition.PCA]
    sample_names: List[str]
    sample_coordinate_df: pd.DataFrame

    def main_workflow(self):
        self.assert_data_structure()
        self.transpose_for_row_features()
        self.set_embedding()
        self.set_sample_names()
        self.fit_transform()

    def assert_data_structure(self):
        assert self.data_structure in self.DATA_STRUCTURES

    def transpose_for_row_features(self):
        if self.data_structure == 'row_features':
            self.df = self.df.transpose()

    def set_sample_names(self):
        self.sample_names = list(self.df.index) \
            if self.data_structure == 'row_features' \
            else list(self.df.columns)

    @abstractmethod
    def set_embedding(self):
        pass

    def fit_transform(self):
        transformed = self.embedding.fit_transform(
            self.df.to_numpy()
        )

        self.sample_coordinate_df = pd.DataFrame(
            data=transformed,
            columns=self.XY_COLUMNS,
            index=self.sample_names)


class PCACore(Core):

    XY_COLUMNS = ('PC 1', 'PC 2')

    embedding: decomposition.PCA
    proportion_explained_series: pd.Series

    def main(
            self,
            df: pd.DataFrame,
            data_structure: str) -> Tuple[pd.DataFrame, pd.Series]:

        self.df = df
        self.data_structure = data_structure

        self.main_workflow()
        self.set_proportion_explained_serise()

        return self.sample_coordinate_df, self.proportion_explained_series

    def set_embedding(self):
        self.embedding = decomposition.PCA(
            n_components=self.N_COMPONENTS,
            copy=True,
            whiten=False,
            svd_solver='auto',
            tol=0.0,
            iterated_power='auto',
            random_state=self.RANDOM_STATE)

    def set_proportion_explained_serise(self):
        self.proportion_explained_series = self.embedding.explained_variance_ratio_


class NMDSCore(Core):

    XY_COLUMNS = ('NMDS 1', 'NMDS 2')
    METRIC = False  # i.e. non-metric MDS
    N_INIT = 10  # number of independent fitting
    EPS = 1e-3  # stress tolerance for convergence

    embedding: manifold.MDS
    stress: float

    def main(
            self,
            df: pd.DataFrame,
            data_structure: str) -> Tuple[pd.DataFrame, float]:

        self.df = df
        self.data_structure = data_structure

        self.main_workflow()
        self.normalize_stress()

        return self.sample_coordinate_df, self.stress

    def set_embedding(self):
        dissimilarity = 'precomputed' \
            if self.data_structure == 'distance_matrix' \
            else 'euclidean'

        self.embedding = manifold.MDS(
            n_components=self.N_COMPONENTS,
            metric=self.METRIC,
            n_init=self.N_INIT,
            max_iter=300,
            verbose=0,
            eps=self.EPS,
            n_jobs=self.threads,
            random_state=self.RANDOM_STATE,
            dissimilarity=dissimilarity)

    def normalize_stress(self):
        """
        Normalize with sum of squared distances
        i.e. Kruskal Stress, or Stress_1
        https://stackoverflow.com/questions/36428205/stress-attribute-sklearn-manifold-mds-python
        """
        dist_mat = self.embedding.dissimilarity_matrix_
        squared_sum = np.sum(dist_mat ** 2) / 2  # diagonal symmetry, thus divide by 2
        self.stress = np.sqrt(self.embedding.stress_ / squared_sum)


class TSNECore(Core):

    XY_COLUMNS = ('t-SNE 1', 't-SNE 2')
    PERPLEXITY = 3.0

    embedding: manifold.TSNE

    def main(
            self,
            df: pd.DataFrame,
            data_structure: str) -> pd.DataFrame:

        self.df = df
        self.data_structure = data_structure

        self.main_workflow()

        return self.sample_coordinate_df

    def set_embedding(self):
        if self.data_structure == 'distance_matrix':
            metric = 'precomputed'
            init = 'random'
        else:
            metric = 'euclidean'
            init = 'pca'

        self.embedding = manifold.TSNE(
            n_components=self.N_COMPONENTS,
            perplexity=self.PERPLEXITY,
            early_exaggeration=12.0,
            learning_rate=200.0,
            n_iter=1000,
            n_iter_without_progress=300,
            min_grad_norm=1e-7,
            metric=metric,
            init=init,
            verbose=1,
            random_state=self.RANDOM_STATE,
            method='barnes_hut',
            angle=0.5,
            n_jobs=self.threads,
            square_distances='legacy'
        )


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
