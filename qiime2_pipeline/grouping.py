import pandas as pd
from .template import Processor


class AddGroupColumn(Processor):

    GROUP_COLUMN: str = 'Group'
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
        assert 'Group' in self.sample_df.columns, \
            f'No "{self.GROUP_COLUMN}" column in {self.sample_sheet}'

    def add_group_column(self):
        self.df = self.df.merge(
            right=self.sample_df[self.GROUP_COLUMN],
            how='left',
            left_index=True,
            right_index=True)
        self.df[self.GROUP_COLUMN] = self.df[self.GROUP_COLUMN].fillna(self.NA_VALUE)

    def reorder_columns(self):
        cols = list(self.df.columns)
        cols = [cols[-1]] + cols[:-1]
        self.df = self.df[cols]
