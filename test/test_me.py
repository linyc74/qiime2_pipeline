from qiime2_pipeline.beta import BetaDiversity
from .setup import TestCase


class MyTest(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test(self):
        BetaDiversity(self.settings).main(
            feature_tabe_qza=f'{self.indir}/dada2-table.qza',
            rooted_tree_qza=f'{self.indir}/rooted-fasttree.qza'
        )
