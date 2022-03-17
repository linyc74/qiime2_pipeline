import pandas as pd
from os.path import exists
from qiime2_pipeline.beta_embedding import PCoA, NMDS, TSNE, ScatterPlot, BatchPCoA, BatchNMDS, BatchTSNE
from .setup import TestCase


class TestPCoA(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        PCoA(self.settings).main(
            distance_matrix_tsv=f'{self.indir}/distance-matrix.tsv',
            group_keywords=['H']
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-proportion-explained.tsv',
        ]:
            self.assertTrue(exists(f))


class TestNMDS(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        NMDS(self.settings).main(
            distance_matrix_tsv=f'{self.indir}/distance-matrix.tsv',
            group_keywords=['H']
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-stress.txt',
        ]:
            self.assertTrue(exists(f))


class TestTSNE(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        TSNE(self.settings).main(
            distance_matrix_tsv=f'{self.indir}/distance-matrix.tsv',
            group_keywords=['H']
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.png',
        ]:
            self.assertTrue(exists(f))


class TestBatchPCoA(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        BatchPCoA(self.settings).main(
            distance_matrix_tsvs=[f'{self.indir}/distance-matrix.tsv'],
            group_keywords=['H']
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-proportion-explained.tsv',
        ]:
            self.assertTrue(exists(f))


class TestBatchNMDS(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        BatchNMDS(self.settings).main(
            distance_matrix_tsvs=[f'{self.indir}/distance-matrix.tsv'],
            group_keywords=['H']
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-stress.txt',
        ]:
            self.assertTrue(exists(f))


class TestBatchTSNE(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        BatchTSNE(self.settings).main(
            distance_matrix_tsvs=[f'{self.indir}/distance-matrix.tsv'],
            group_keywords=['H']
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.png',
        ]:
            self.assertTrue(exists(f))


class TestScatterPlot(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        sample_coordinate_df = pd.read_csv(
            f'{self.indir}/sample-coordinate.tsv',
            index_col=0,
            sep='\t'
        )
        ScatterPlot(self.settings).main(
            sample_coordinate_df=sample_coordinate_df,
            x_column='PC1',
            y_column='PC2',
            hue_column='Group',
            output_png=f'{self.outdir}/sample-coordinate.png'
        )
