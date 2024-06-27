import pandas as pd
from qiime2_pipeline.sample_sheet import TranscribeSampleSheet
from .setup import TestCase


class TestTranscribeSampleSheet(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_csv(self):
        actual = TranscribeSampleSheet(self.settings).main(sample_sheet=f'{self.indir}/sample-sheet.csv')
        self.assertDataFrameEqual(
            pd.read_csv(f'{self.indir}/expected.csv'),
            pd.read_csv(actual)
        )

    def test_tsv(self):
        actual = TranscribeSampleSheet(self.settings).main(sample_sheet=f'{self.indir}/sample-sheet.tsv')
        self.assertDataFrameEqual(
            pd.read_csv(f'{self.indir}/expected.csv'),
            pd.read_csv(actual)
        )

    def test_xlsx(self):
        actual = TranscribeSampleSheet(self.settings).main(sample_sheet=f'{self.indir}/sample-sheet.xlsx')
        self.assertDataFrameEqual(
            pd.read_csv(f'{self.indir}/expected.csv'),
            pd.read_csv(actual)
        )

    def test_remove_empty_rows_and_columns(self):
        actual = TranscribeSampleSheet(self.settings).main(sample_sheet=f'{self.indir}/faulty-sample-sheet.csv')
        self.assertDataFrameEqual(
            pd.read_csv(f'{self.indir}/expected.csv'),
            pd.read_csv(actual)
        )
