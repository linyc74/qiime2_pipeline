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
            feature_table_tsv=f'{self.indir}/feature-table.tsv',
            rooted_tree_qza=f'{self.indir}/fasttree-rooted.qza',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0),  'green'],
        )


#


class TestRunAllBetaMetricsToTsvs(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = RunAllBetaMetricsToDistanceMatrixTsvs(self.settings).main(
            feature_table_qza=f'{self.indir}/feature-table.qza',
            rooted_tree_qza=f'{self.indir}/fasttree-rooted.qza'
        )
        expected = [
            f'{self.outdir}/beta-diversity/jaccard.tsv',
            f'{self.outdir}/beta-diversity/euclidean.tsv',
            f'{self.outdir}/beta-diversity/braycurtis.tsv',
            f'{self.outdir}/beta-diversity/weighted_unifrac.tsv',
            f'{self.outdir}/beta-diversity/weighted_normalized_unifrac.tsv',
            f'{self.outdir}/beta-diversity/generalized_unifrac.tsv',
            f'{self.outdir}/beta-diversity/unweighted_unifrac.tsv',
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
        name = 'jaccard'
        PCoAProcess(self.settings).main(
            tsv=f'{self.indir}/distance-matrix-tsvs/{name}.tsv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0),  'green'],
        )
        for f in [
            f'{self.outdir}/beta-embedding/{name}-pcoa-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/{name}-pcoa-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/{name}-pcoa-sample-coordinate.pdf',
            f'{self.outdir}/beta-embedding/{name}-pcoa-proportion-explained.tsv',
        ]:
            self.assertTrue(exists(f))


class TestTSNEProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        name = 'jaccard'
        TSNEProcess(self.settings).main(
            tsv=f'{self.indir}/distance-matrix-tsvs/{name}.tsv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)],
        )
        for f in [
            f'{self.outdir}/beta-embedding/{name}-tsne-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/{name}-tsne-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/{name}-tsne-sample-coordinate.pdf',
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
            distance_matrix_tsvs=[
                f'{self.indir}/distance-matrix-tsvs/braycurtis.tsv',
                f'{self.indir}/distance-matrix-tsvs/euclidean.tsv',
                f'{self.indir}/distance-matrix-tsvs/jaccard.tsv',
                f'{self.indir}/distance-matrix-tsvs/unweighted_unifrac.tsv',
                f'{self.indir}/distance-matrix-tsvs/weighted_unifrac.tsv',
                f'{self.indir}/distance-matrix-tsvs/weighted_normalized_unifrac.tsv',
                f'{self.indir}/distance-matrix-tsvs/generalized_unifrac.tsv',
            ],
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0),  'green'],
        )


class TestBatchTSNEProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        BatchTSNEProcess(self.settings).main(
            distance_matrix_tsvs=[
                f'{self.indir}/distance-matrix-tsvs/braycurtis.tsv',
                f'{self.indir}/distance-matrix-tsvs/euclidean.tsv',
                f'{self.indir}/distance-matrix-tsvs/jaccard.tsv',
                f'{self.indir}/distance-matrix-tsvs/unweighted_unifrac.tsv',
                f'{self.indir}/distance-matrix-tsvs/weighted_unifrac.tsv',
                f'{self.indir}/distance-matrix-tsvs/weighted_normalized_unifrac.tsv',
                f'{self.indir}/distance-matrix-tsvs/generalized_unifrac.tsv',
            ],
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0),  'green'],
        )
