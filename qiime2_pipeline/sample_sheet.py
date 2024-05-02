import pandas as pd
from .template import Processor


class TranscribeSampleSheet(Processor):

    sample_sheet: str
    df: pd.DataFrame
    output_csv: str

    def main(self, sample_sheet: str) -> str:

        self.sample_sheet = sample_sheet

        self.output_csv = f'{self.workdir}/sample-sheet.csv'

        if self.sample_sheet.endswith('.csv'):
            self.df = pd.read_csv(self.sample_sheet)

        elif self.sample_sheet.endswith('.xlsx'):
            self.df = pd.read_excel(self.sample_sheet)

        else:  # default to tab-delimited for .txt or .tsv
            self.df = pd.read_csv(self.sample_sheet, sep='\t')

        self.df.to_csv(self.output_csv, index=False)

        return self.output_csv
