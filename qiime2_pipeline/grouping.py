import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from .template import Processor


GROUP_COLUMN = 'Group'


class AddGroupColumn(Processor):

    NA_VALUE: str = 'None'  # Don't use 'NA', which would make dtype = `float` but not `str`, tricky for testing

    df: pd.DataFrame
    sample_sheet: str

    sample_df: pd.DataFrame

    def main(
            self,
            df: pd.DataFrame,
            sample_sheet: str) -> pd.DataFrame:

        self.df = df.copy(deep=True)
        self.sample_sheet = sample_sheet

        self.read_sample_sheet()
        self.add_group_column()
        self.reorder_columns()

        return self.df

    def read_sample_sheet(self):
        self.sample_df = pd.read_csv(self.sample_sheet, index_col=0)
        assert GROUP_COLUMN in self.sample_df.columns, \
            f'No "{GROUP_COLUMN}" column in {self.sample_sheet}'

    def add_group_column(self):
        self.df = self.df.merge(
            right=self.sample_df[GROUP_COLUMN],
            how='left',
            left_index=True,
            right_index=True)
        self.df[GROUP_COLUMN] = self.df[GROUP_COLUMN].fillna(self.NA_VALUE)

    def reorder_columns(self):
        cols = list(self.df.columns)
        cols = [cols[-1]] + cols[:-1]
        self.df = self.df[cols]


class TagGroupNamesOnSampleColumns(Processor):

    df: pd.DataFrame
    sample_sheet: str

    sample_df: pd.DataFrame

    def main(
            self,
            df: pd.DataFrame,
            sample_sheet: str) -> pd.DataFrame:

        self.df = df.copy(deep=True)
        self.sample_sheet = sample_sheet

        self.read_sample_sheet()
        self.rename_columns()

        return self.df

    def read_sample_sheet(self):
        self.sample_df = pd.read_csv(self.sample_sheet, index_col=0)
        assert GROUP_COLUMN in self.sample_df.columns, \
            f'No "{GROUP_COLUMN}" column in {self.sample_sheet}'

    def rename_columns(self):
        for col in self.df.columns:
            if col in self.sample_df.index:
                group = self.sample_df.loc[col, GROUP_COLUMN]
                self.df = self.df.rename(columns={col: f'[{group}] {col}'})


class GetColors(Processor):

    sample_sheet: str
    colormap: str
    invert_colors: bool

    def main(
            self,
            sample_sheet: str,
            colormap: str,
            invert_colors: bool) -> list:

        self.sample_sheet = sample_sheet
        self.colormap = colormap
        self.invert_colors = invert_colors

        df = pd.read_csv(self.sample_sheet, index_col=0)
        n_groups = len(df[GROUP_COLUMN].unique())

        if ',' in self.colormap:
            names = self.colormap.split(',')
            if len(names) != n_groups:
                self.logger.info(f'WARNING! Number of colors "{self.colormap}" does not match number of groups ({n_groups})')
            colors = [to_rgba(n) for n in names]
        else:
            cmap = plt.colormaps[self.colormap]
            colors = [cmap(i) for i in range(n_groups)]

        if self.invert_colors:
            colors = colors[::-1]

        return colors
