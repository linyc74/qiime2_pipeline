import pandas as pd
from typing import List
from .template import Processor


class AddGroupColumn(Processor):

    GROUP_COLUMN: str = 'Group'
    NA_VALUE: str = 'None'  # Don't use 'NA', which would make dtype = `float` but not `str`, tricky for testing

    df: pd.DataFrame
    group_keywords: List[str]

    def main(
            self,
            df: pd.DataFrame,
            group_keywords: List[str]) -> pd.DataFrame:

        self.df = df.copy(deep=True)
        self.group_keywords = group_keywords

        self.add_group_column_using_index()
        self.reorder_columns()

        return self.df

    def add_group_column_using_index(self):
        for idx in self.df.index:
            self.df.loc[idx, self.GROUP_COLUMN] = self.__to_group(idx=idx)

    def __to_group(self, idx: str) -> str:
        for k in self.group_keywords:
            if k in idx:
                return k
        return self.NA_VALUE

    def reorder_columns(self):
        cols = list(self.df.columns)
        cols = [cols[-1]] + cols[:-1]
        self.df = self.df[cols]
