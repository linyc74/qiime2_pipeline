import os
import pandas as pd
import seaborn as sns
from sklearn import manifold, decomposition
from typing import Tuple, Union
from matplotlib.axes import Axes
from matplotlib import pyplot as plt
from abc import ABC, abstractmethod
from .tools import edit_fpath
from .template import Processor
from .grouping import AddGroupColumn


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
    sample_coordinate_df: pd.DataFrame

    def main_workflow(self):
        self.assert_data_structure()
        self.transpose_for_row_features()
        self.set_embedding()
        self.fit_transform()

    def assert_data_structure(self):
        assert self.data_structure in self.DATA_STRUCTURES

    def transpose_for_row_features(self):
        if self.data_structure == 'row_features':
            self.df = self.df.transpose()

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
            index=list(self.df.index)
        )


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
        self.proportion_explained_series = pd.Series(self.embedding.explained_variance_ratio_)


class TSNECore(Core):

    XY_COLUMNS = ('t-SNE 1', 't-SNE 2')
    PERPLEXITY = 5.0

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


class EmbeddingProcessTemplate(Processor, ABC):

    NAME: str
    XY_COLUMNS: Tuple[str, str]
    DSTDIR_NAME: str
    GROUP_COLUMN: str = AddGroupColumn.GROUP_COLUMN

    tsv: str
    sample_sheet: str
    colormap: str

    df: pd.DataFrame
    sample_coordinate_df: pd.DataFrame
    dstdir: str

    @abstractmethod
    def main(
            self,
            tsv: str,
            sample_sheet: str,
            colormap: str):
        pass

    def run_main_workflow(self):
        self.read_tsv()

        self.preprocessing()
        self.embedding()

        self.add_group_column()
        self.make_dstdir()
        self.write_sample_coordinate()
        self.plot_sample_coordinate()

    def read_tsv(self):
        self.df = pd.read_csv(
            self.tsv,
            sep='\t',
            index_col=0)

    @abstractmethod
    def preprocessing(self):
        pass

    @abstractmethod
    def embedding(self):
        pass

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def add_group_column(self):
        self.sample_coordinate_df = AddGroupColumn(self.settings).main(
            df=self.sample_coordinate_df,
            sample_sheet=self.sample_sheet)

    def write_sample_coordinate(self):
        tsv = self.__get_sample_coordinate_fpath(suffix='.tsv')
        self.sample_coordinate_df.to_csv(tsv, sep='\t')

    def plot_sample_coordinate(self):
        output_prefix = self.__get_sample_coordinate_fpath(suffix='')
        ScatterPlot(self.settings).main(
            sample_coordinate_df=self.sample_coordinate_df,
            x_column=self.XY_COLUMNS[0],
            y_column=self.XY_COLUMNS[1],
            hue_column=self.GROUP_COLUMN,
            colormap=self.colormap,
            output_prefix=output_prefix)

    def __get_sample_coordinate_fpath(self, suffix: str) -> str:
        name = self.NAME.lower().replace('-', '')
        return edit_fpath(
            fpath=self.tsv,
            old_suffix='.tsv',
            new_suffix=f'-{name}-sample-coordinate{suffix}',
            dstdir=self.dstdir
        )


class ScatterPlot(Processor):

    FIGSIZE = (8, 8)
    DPI = 600

    sample_coordinate_df: pd.DataFrame
    x_column: str
    y_column: str
    group_column: str
    colormap: str
    output_prefix: str

    ax: Axes

    def main(
            self,
            sample_coordinate_df: pd.DataFrame,
            x_column: str,
            y_column: str,
            hue_column: str,
            colormap: str,
            output_prefix: str):

        self.sample_coordinate_df = sample_coordinate_df
        self.x_column = x_column
        self.y_column = y_column
        self.group_column = hue_column
        self.colormap = colormap
        self.output_prefix = output_prefix

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
            hue=self.group_column,
            palette=self.colormap)

    def label_points(self):
        df = self.sample_coordinate_df
        for sample_name in df.index:
            self.ax.text(
                x=df.loc[sample_name, self.x_column],
                y=df.loc[sample_name, self.y_column],
                s=sample_name
            )

    def save_figure(self):
        for ext in ['pdf', 'png']:
            plt.savefig(f'{self.output_prefix}.{ext}', dpi=self.DPI)
        plt.close()
