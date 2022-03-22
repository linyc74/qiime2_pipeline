import pandas as pd
from .setup import TestCase
from qiime2_pipeline.normalization import CountNormalization


class TestCountNormalization(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        df = pd.read_csv(f'{self.indir}/count-table.tsv', sep='\t', index_col=0)
        actual = CountNormalization(self.settings).main(
            df=df,
            log_pseudocount=True,
            by_sample_reads=True,
            sample_reads_unit=10000
        )
        expected = pd.read_csv(f'{self.indir}/count-table-normalized.tsv', sep='\t', index_col=0)
        self.assertDataFrameEqual(expected, actual)
