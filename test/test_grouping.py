import pandas as pd
from .setup import TestCase
from qiime2_pipeline.grouping import Grouping


class TestGrouping(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        outdf = Grouping(self.settings).main(
            indf=pd.read_csv(f'{self.indir}/indf.csv', index_col=0),
            group_keywords=['H', 'O']
        )
        self.assertDataFrameEqual(
            first=pd.read_csv(f'{self.indir}/outdf.csv', index_col=0),
            second=outdf
        )
