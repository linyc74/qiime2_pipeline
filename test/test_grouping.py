import pandas as pd
from .setup import TestCase
from qiime2_pipeline.grouping import AddGroupColumn, TagGroupNamesOnSampleColumns, GetColors


class TestAddGroupColumn(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = AddGroupColumn(self.settings).main(
            df=pd.read_csv(f'{self.indir}/indf.csv', index_col=0),
            sample_sheet=f'{self.indir}/sample-sheet.csv'
        )
        expected = pd.read_csv(f'{self.indir}/outdf.csv', index_col=0)
        self.assertDataFrameEqual(expected, actual)


class TestTagGroupNameOnSampleColumns(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = TagGroupNamesOnSampleColumns(self.settings).main(
            df=pd.read_csv(f'{self.indir}/phylum-table.tsv', index_col=0, sep='\t'),
            sample_sheet=f'{self.indir}/sample-sheet.csv'
        )
        expected = pd.read_csv(f'{self.indir}/phylum-table-tagged.tsv', index_col=0, sep='\t')
        self.assertDataFrameEqual(expected, actual)


class TestGetColors(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = GetColors(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colormap='Set1',
            invert_colors=True
        )
        expected = [  # two colors
            (0.21568627450980393, 0.49411764705882355, 0.7215686274509804, 1.0),
            (0.8941176470588236, 0.10196078431372549, 0.10980392156862745, 1.0),
        ]
        self.assertListEqual(expected, actual)
