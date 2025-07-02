import os
import pandas as pd
import seaborn as sns
from sklearn import manifold, decomposition
from matplotlib.axes import Axes
from matplotlib import pyplot as plt
from abc import ABC, abstractmethod
from typing import Tuple, Union, List
from .utils import edit_fpath
from .template import Processor
from .grouping import GROUP_COLUMN, AddGroupColumn


class Core(Processor, ABC):

    XY_COLUMNS: List[str]
    DATA_STRUCTURES = [
        'distance_matrix',
        'row_features',
        'column_features',
    ]
    N_COMPONENTS = 2
    RANDOM_STATE = 1  # to ensure reproducible result

    df: pd.DataFrame
    data_structure: str

    embedding: Union[decomposition.PCA]
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

    XY_COLUMNS = ['PC 1', 'PC 2']

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


class EmbeddingProcessTemplate(Processor, ABC):

    NAME: str
    XY_COLUMNS: Tuple[str, str]
    DSTDIR_NAME: str
    GROUP_COLUMN: str = GROUP_COLUMN

    tsv: str
    sample_sheet: str
    colors: list

    df: pd.DataFrame
    sample_coordinate_df: pd.DataFrame
    dstdir: str

    @abstractmethod
    def main(
            self,
            tsv: str,
            sample_sheet: str,
            colors: list):
        pass

    def run_main_workflow(self):
        self.read_tsv()

        self.preprocessing()
        self.embedding()

        self.reorder_sample_coordinate_df()
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

    def reorder_sample_coordinate_df(self):
        # reorder sample ids according to sample sheet, to keep the color order consistent
        # for sample ids not in sample sheet, append them at the end
        all_ids = self.sample_coordinate_df.index.tolist()
        sample_sheet_ids = pd.read_csv(self.sample_sheet, index_col=0).index.tolist()
        other_ids = [s for s in all_ids if s not in sample_sheet_ids]
        self.sample_coordinate_df = self.sample_coordinate_df.loc[sample_sheet_ids + other_ids]

    def add_group_column(self):
        self.sample_coordinate_df = AddGroupColumn(self.settings).main(
            df=self.sample_coordinate_df,
            sample_sheet=self.sample_sheet)

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

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
            colors=self.colors,
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

    sample_coordinate_df: pd.DataFrame
    x_column: str
    y_column: str
    group_column: str
    colors: list
    output_prefix: str

    figsize: Tuple[float, float]
    point_size: float
    marker_edge_color: str
    line_width: float
    fontsize: int
    dpi: int

    ax: Axes

    def main(
            self,
            sample_coordinate_df: pd.DataFrame,
            x_column: str,
            y_column: str,
            hue_column: str,
            colors: list,
            output_prefix: str):

        self.sample_coordinate_df = sample_coordinate_df.copy()
        self.x_column = x_column
        self.y_column = y_column
        self.group_column = hue_column
        self.colors = colors
        self.output_prefix = output_prefix

        self.set_figsize()
        self.set_parameters()
        self.init_figure()
        self.scatterplot()
        self.label_points()
        self.save_figure()

    def set_figsize(self):
        df, group = self.sample_coordinate_df, self.group_column
        max_legend_chrs = get_max_str_length(df[group])
        n_groups = len(df[group].unique())

        self.figsize = GetFigsize(self.settings).main(
            max_legend_chrs=max_legend_chrs, n_groups=n_groups)

    def set_parameters(self):
        if self.settings.for_publication:
            self.point_size = 20.
            self.marker_edge_color = 'white'
            self.line_width = 0.5
            self.fontsize = 7
            self.dpi = 600
        else:
            self.point_size = 30.
            self.marker_edge_color = 'white'
            self.line_width = 1.0
            self.fontsize = 10
            self.dpi = 300

    def init_figure(self):
        plt.rcParams['font.size'] = self.fontsize
        plt.rcParams['axes.linewidth'] = self.line_width
        plt.figure(figsize=self.figsize, dpi=self.dpi)

    def scatterplot(self):
        self.ax = sns.scatterplot(
            data=self.sample_coordinate_df,
            x=self.x_column,
            y=self.y_column,
            s=self.point_size,
            edgecolor='white',
            hue=self.group_column,
            palette=self.colors
        )
        plt.gca().xaxis.set_tick_params(width=self.line_width)
        plt.gca().yaxis.set_tick_params(width=self.line_width)
        legend = plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        legend.set_frame_on(False)

    def label_points(self):
        if self.settings.for_publication:
            return
        df = self.sample_coordinate_df
        for sample_name in df.index:
            self.ax.text(
                x=df.loc[sample_name, self.x_column],
                y=df.loc[sample_name, self.y_column],
                s=sample_name,
                alpha=0.25,
                fontsize=6
            )

    def save_figure(self):
        plt.tight_layout()
        for ext in ['pdf', 'png']:
            plt.savefig(f'{self.output_prefix}.{ext}', dpi=self.dpi)
        plt.close()


def get_max_str_length(series: pd.Series) -> int:
    return series.astype(str).apply(len).max()


class GetFigsize(Processor):

    def main(self, max_legend_chrs: int, n_groups: int) -> Tuple[float, float]:

        if self.settings.for_publication:
            base_width = 7.5 / 2.54
            chr_width = 0.15 / 2.54
            base_height = 6 / 2.54
            line_height = 0.4 / 2.54
        else:
            base_width = 14 / 2.54
            chr_width = 0.218 / 2.54
            base_height = 12 / 2.54
            line_height = 0.5 / 2.54

        w = base_width + (max_legend_chrs * chr_width)
        h = max(base_height, n_groups * line_height)

        return w, h
