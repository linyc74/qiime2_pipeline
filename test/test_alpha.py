import pandas as pd
from os.path import exists
from .setup import TestCase
from qiime2_pipeline.alpha import AlphaDiversity, Plot, Statistics


class TestAlphaDiversity(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_unpaired(self):
        self.settings.for_publication = True
        AlphaDiversity(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            feature_table_qza=f'{self.indir}/feature-table.qza',
            subject_column=None,
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
            f'{self.outdir}/alpha-diversity/observed_features.png',
            f'{self.outdir}/alpha-diversity/chao1.png',
            f'{self.outdir}/alpha-diversity/shannon_entropy.png',
            f'{self.outdir}/alpha-diversity/simpson.png',
        ]:
            self.assertTrue(exists(file))

    def test_paired(self):
        self.settings.for_publication = True
        AlphaDiversity(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            feature_table_qza=f'{self.indir}/feature-table.qza',
            subject_column='Subject',
            alpha_metrics=[
                'observed_features',
            ],
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )
        for file in [
            f'{self.outdir}/alpha-diversity/alpha-diversity.csv',
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
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)],
            subject_column=None
        )

    def test_for_publication(self):
        self.settings.for_publication = True
        Plot(self.settings).main(
            df=pd.read_csv(f'{self.indir}/alpha-diversity-unpaired.csv', index_col=0),
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)],
            subject_column=None
        )
    
    def test_with_subject_column(self):
        self.settings.for_publication = True
        Plot(self.settings).main(
            df=pd.read_csv(f'{self.indir}/alpha-diversity-paired.csv', index_col=0),
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)],
            subject_column='Subject'
        )


class TestStatistics(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_unpaired(self):
        Statistics(self.settings).main(
            df=pd.read_csv(f'{self.indir}/alpha-diversity-unpaired.csv', index_col=0),
            subject_column=None)

    def test_paired(self):
        Statistics(self.settings).main(
            df=pd.read_csv(f'{self.indir}/alpha-diversity-paired.csv', index_col=0),
            subject_column='Subject')
