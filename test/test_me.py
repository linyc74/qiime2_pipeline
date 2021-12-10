from .setup import TestCase
from qiime2_pipeline.beta import BetaDiversity


class MyTest(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_me(self):
        BetaDiversity(self.settings).main(
            feature_table_qza=f'{self.indir}/dada2-table.qza',
            rooted_tree_qza=f'{self.indir}/rooted-fasttree.qza',
        )
