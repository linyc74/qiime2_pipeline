from os.path import exists
from .setup import TestCase
from qiime2_pipeline.beta_qiime import QiimeBetaDiversity, RunAllBetaMetricsToDistanceMatrixTsvs, \
    PCoAProcess, TSNEProcess, BatchPCoAProcess, BatchTSNEProcess


class TestQiimeBetaDiversity(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        QiimeBetaDiversity(self.settings).main(
            feature_table_qza=f'{self.indir}/labeled-feature-table-normalized.qza',
            rooted_tree_qza=f'{self.indir}/mafft-aligned-sequences-masked.qza',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colormap='Set1',
            invert_colors=True,
        )


#


class TestRunAllBetaMetricsToTsvs(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = RunAllBetaMetricsToDistanceMatrixTsvs(self.settings).main(
            feature_table_qza=f'{self.indir}/labeled-feature-table-normalized.qza',
            rooted_tree_qza=f'{self.indir}/mafft-aligned-sequences-masked.qza'
        )
        expected = [
            f'{self.outdir}/beta-diversity/jaccard.tsv',
            f'{self.outdir}/beta-diversity/braycurtis.tsv',
        ]
        for e, a in zip(expected, actual):
            self.assertFileExists(expected=e, actual=a)


#


class TestPCoAProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        PCoAProcess(self.settings).main(
            tsv=f'{self.indir}/distance-matrix.tsv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colormap='Set1',
            invert_colors=True,
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.pdf',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-proportion-explained.tsv',
        ]:
            self.assertTrue(exists(f))


class TestTSNEProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        TSNEProcess(self.settings).main(
            tsv=f'{self.indir}/distance-matrix.tsv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colormap='Set1',
            invert_colors=True,
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.pdf',
        ]:
            self.assertTrue(exists(f))


#


class TestBatchPCoAProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        BatchPCoAProcess(self.settings).main(
            distance_matrix_tsvs=[f'{self.indir}/distance-matrix.tsv'],
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colormap='Set1',
            invert_colors=True,
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.pdf',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-proportion-explained.tsv',
        ]:
            self.assertTrue(exists(f))


class TestBatchTSNEProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        BatchTSNEProcess(self.settings).main(
            distance_matrix_tsvs=[f'{self.indir}/distance-matrix.tsv'],
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colormap='Set1',
            invert_colors=True,
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.pdf',
        ]:
            self.assertTrue(exists(f))
