from os.path import exists
from .setup import TestCase
from qiime2_pipeline.beta_qiime import QiimeBetaDiversity, RunAllBetaMetricsToDistanceMatrixTsvs, \
    PCoAProcess, NMDSProcess, TSNEProcess, BatchPCoAProcess, BatchNMDSProcess, BatchTSNEProcess


class TestQiimeBetaDiversity(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        QiimeBetaDiversity(self.settings).main(
            feature_table_qza=f'{self.indir}/labeled-feature-table-normalized.qza',
            rooted_tree_qza=f'{self.indir}/mafft-aligned-sequences-masked.qza',
            group_keywords=['H']
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
            f'{self.outdir}/beta-diversity/cosine.tsv',
            f'{self.outdir}/beta-diversity/correlation.tsv',
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
            group_keywords=['H']
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.pdf',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-proportion-explained.tsv',
        ]:
            self.assertTrue(exists(f))


class TestNMDSProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        NMDSProcess(self.settings).main(
            tsv=f'{self.indir}/distance-matrix.tsv',
            group_keywords=['H']
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-sample-coordinate.pdf',
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-stress.txt',
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
            group_keywords=['H']
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
            group_keywords=['H']
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-sample-coordinate.pdf',
            f'{self.outdir}/beta-embedding/distance-matrix-pcoa-proportion-explained.tsv',
        ]:
            self.assertTrue(exists(f))


class TestBatchNMDSProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        BatchNMDSProcess(self.settings).main(
            distance_matrix_tsvs=[f'{self.indir}/distance-matrix.tsv'],
            group_keywords=['H']
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-sample-coordinate.pdf',
            f'{self.outdir}/beta-embedding/distance-matrix-nmds-stress.txt',
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
            group_keywords=['H']
        )
        for f in [
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/distance-matrix-tsne-sample-coordinate.pdf',
        ]:
            self.assertTrue(exists(f))
