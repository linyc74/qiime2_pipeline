import pandas as pd
from os.path import exists
from .setup import TestCase
from qiime2_pipeline.alpha import AlphaDiversity, Plot


class TestAlphaDiversity(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        self.settings.for_publication = True
        AlphaDiversity(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            feature_table_qza=f'{self.indir}/feature-table.qza',
            alpha_metrics=[
                'observed_features',
                'chao1',
                'shannon',
                'simpson'
            ],
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )
        for file in [
            f'{self.outdir}/alpha-diversity/alpha-diversity.csv',
            f'{self.outdir}/alpha-diversity/chao1.png',
            f'{self.outdir}/alpha-diversity/shannon_entropy.png',
            f'{self.outdir}/alpha-diversity/simpson.png',
            f'{self.outdir}/alpha-diversity/observed_features.png',
        ]:
            self.assertTrue(exists(file))


class TestPlot(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_not_for_publication(self):
        self.settings.for_publication = False
        Plot(self.settings).main(
            df=pd.read_csv(f'{self.indir}/alpha-diversity-346-samples.csv', index_col=0),
            dstdir=self.outdir,
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )

    def test_for_publication(self):
        self.settings.for_publication = True
        Plot(self.settings).main(
            df=pd.read_csv(f'{self.indir}/alpha-diversity.csv', index_col=0),
            dstdir=self.outdir,
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )
