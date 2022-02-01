import pandas as pd
from typing import List
from .template import Processor, Settings


class Grouping(Processor):

    GROUP = 'Group'
    NA_VALUE = 'None'  # Don't use 'NA', which would make dtype = float

    indf: pd.DataFrame
    group_keywords: List[str]

    outdf: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            indf: pd.DataFrame,
            group_keywords: List[str]) -> pd.DataFrame:

        self.indf = indf
        self.group_keywords = group_keywords

        self.outdf = self.indf.copy()
        self.add_group_column()
        self.reorder_columns()

        return self.outdf

    def add_group_column(self):
        for idx in self.outdf.index:
            self.outdf.loc[idx, self.GROUP] = self.__idx_to_group(idx=idx)

    def __idx_to_group(self, idx: str) -> str:
        for k in self.group_keywords:
            if k in idx:
                return k
        return self.NA_VALUE

    def reorder_columns(self):
        cols = list(self.outdf.columns)
        cols = [cols[-1]] + cols[:-1]
        self.outdf = self.outdf[cols]
