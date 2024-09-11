from os.path import exists
from .setup import TestCase
from qiime2_pipeline.alpha_rarefaction import AlphaRarefaction


class TestAlphaRarefaction(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        AlphaRarefaction(self.settings).main(
            feature_table_qza=f'{self.indir}/feature-table.qza',
        )
        self.assertTrue(exists(f'{self.outdir}/rarefaction.qzv'))
