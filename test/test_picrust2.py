from qiime2_pipeline.picrust2 import PICRUSt2
from .setup import TestCase


class TestPICRUSt2(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test_main(self):
        actual = PICRUSt2(self.settings).main(
            labeled_feature_sequence_fa=f'{self.indir}/labeled-feature-sequence.fa',
            labeled_feature_table_tsv=f'{self.indir}/labeled-feature-table.tsv'
        )
        expected = f'{self.outdir}/picrust2-pathway-table.tsv'
        self.assertFileEqual(expected, actual)
