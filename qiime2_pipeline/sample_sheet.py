import numpy as np
import pandas as pd
from .template import Processor


class TranscribeSampleSheet(Processor):

    sample_sheet: str
    df: pd.DataFrame
    output_csv: str

    def main(self, sample_sheet: str) -> str:

        self.sample_sheet = sample_sheet

        self.read_file()
        self.clean_up_empty_rows_and_columns()
        self.save_csv()

        return self.output_csv

    def read_file(self):
        if self.sample_sheet.endswith('.csv'):
            self.df = pd.read_csv(self.sample_sheet)
        elif self.sample_sheet.endswith('.xlsx'):
            self.df = pd.read_excel(self.sample_sheet)
        else:  # default to tab-delimited for .txt or .tsv
            self.df = pd.read_csv(self.sample_sheet, sep='\t')

    def clean_up_empty_rows_and_columns(self):
        for column in self.df.columns:
            if is_empty(self.df[column]):
                self.df.drop(columns=column, inplace=True)
        for idx in self.df.index:
            if is_empty(self.df.loc[idx]):
                self.df.drop(index=idx, inplace=True)

    def save_csv(self):
        self.output_csv = f'{self.workdir}/sample-sheet.csv'
        self.df.to_csv(self.output_csv, index=False)


def is_empty(series: pd.Series) -> bool:
    for value in series:
        if type(value) is str:
            if value.strip(' ') == '':
                value = np.nan  # convert empty string of any length to NaN
        if pd.notna(value):
            return False  # found something not NaN --> not empty
    return True
