from os.path import exists
from .setup import TestCase
from qiime2_pipeline.alpha import AlphaDiversity


class TestAlphaDiversity(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        AlphaDiversity(self.settings).main(
            feature_table_qza=f'{self.indir}/feature-table.qza',
            group_keywords=['H'],
            alpha_metrics=[]  # empty list -> all alpha metrics
        )
        for file in [
            f'{self.outdir}/alpha-diversity/alpha-diversity.csv',
            f'{self.outdir}/alpha-diversity/chao1.png',
            f'{self.outdir}/alpha-diversity/shannon_entropy.png',
            f'{self.outdir}/alpha-diversity/gini_index.png',
            f'{self.outdir}/alpha-diversity/mcintosh_e.png',
            f'{self.outdir}/alpha-diversity/pielou_evenness.png',
            f'{self.outdir}/alpha-diversity/simpson.png',
            f'{self.outdir}/alpha-diversity/observed_features.png',
            f'{self.outdir}/alpha-diversity/fisher_alpha.png',
        ]:
            self.assertTrue(exists(file))
