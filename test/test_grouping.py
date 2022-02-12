import pandas as pd
from .setup import TestCase
from qiime2_pipeline.grouping import AddGroupColumn


class TestAddGroupColumn(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = AddGroupColumn(self.settings).main(
            df=pd.read_csv(f'{self.indir}/indf.csv', index_col=0),
            group_keywords=['H', 'O']
        )
        expected = pd.read_csv(f'{self.indir}/outdf.csv', index_col=0)
        self.assertDataFrameEqual(expected, actual)
