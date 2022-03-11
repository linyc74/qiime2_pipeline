from .setup import TestCase
from qiime2_pipeline.beta import BetaDiversity


class TestBetaDiversity(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        BetaDiversity(self.settings).main(
            feature_table_qza=f'{self.indir}/labeled-feature-table-normalized.qza',
            rooted_tree_qza=f'{self.indir}/mafft-aligned-sequences-masked.qza'
        )
