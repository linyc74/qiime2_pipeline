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
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            feature_table_qza=f'{self.indir}/feature-table.qza',
            alpha_metrics=[
                'observed_features',
                'chao1',
                'shannon',
                'simpson'
            ],
            colormap='Set1',
            invert_colors=True,
        )
        for file in [
            f'{self.outdir}/alpha-diversity/alpha-diversity.csv',
            f'{self.outdir}/alpha-diversity/chao1.png',
            f'{self.outdir}/alpha-diversity/shannon_entropy.png',
            f'{self.outdir}/alpha-diversity/simpson.png',
            f'{self.outdir}/alpha-diversity/observed_features.png',
        ]:
            self.assertTrue(exists(file))
