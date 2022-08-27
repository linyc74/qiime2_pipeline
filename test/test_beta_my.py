from os.path import exists
from qiime2_pipeline.beta_my import MyBetaDiversity, PCAProcess, NMDSProcess, TSNEProcess
from .setup import TestCase


class TestMyBetaDiversity(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        MyBetaDiversity(self.settings).main(
            feature_table_tsv=f'{self.indir}/feature-table.tsv',
            group_keywords=['H', 'O']
        )
        for f in [
            f'{self.outdir}/feature-embedding/feature-table-pca-sample-coordinate.tsv',
            f'{self.outdir}/feature-embedding/feature-table-pca-sample-coordinate.png',
            f'{self.outdir}/feature-embedding/feature-table-pca-sample-coordinate.pdf',
            f'{self.outdir}/feature-embedding/feature-table-pca-proportion-explained.tsv',

            f'{self.outdir}/feature-embedding/feature-table-nmds-sample-coordinate.tsv',
            f'{self.outdir}/feature-embedding/feature-table-nmds-sample-coordinate.png',
            f'{self.outdir}/feature-embedding/feature-table-nmds-sample-coordinate.pdf',
            f'{self.outdir}/feature-embedding/feature-table-nmds-stress.txt',

            f'{self.outdir}/feature-embedding/feature-table-tsne-sample-coordinate.tsv',
            f'{self.outdir}/feature-embedding/feature-table-tsne-sample-coordinate.png',
            f'{self.outdir}/feature-embedding/feature-table-tsne-sample-coordinate.pdf',
        ]:
            self.assertTrue(exists(f))


class TestPCAProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        PCAProcess(self.settings).main(
            tsv=f'{self.indir}/feature-table.tsv',
            group_keywords=['H', 'O']
        )
        for f in [
            f'{self.outdir}/feature-embedding/feature-table-pca-sample-coordinate.tsv',
            f'{self.outdir}/feature-embedding/feature-table-pca-sample-coordinate.png',
            f'{self.outdir}/feature-embedding/feature-table-pca-sample-coordinate.pdf',
            f'{self.outdir}/feature-embedding/feature-table-pca-proportion-explained.tsv',
        ]:
            self.assertTrue(exists(f))


class TestNMDSProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        NMDSProcess(self.settings).main(
            tsv=f'{self.indir}/feature-table.tsv',
            group_keywords=['H', 'O']
        )
        for f in [
            f'{self.outdir}/feature-embedding/feature-table-nmds-sample-coordinate.tsv',
            f'{self.outdir}/feature-embedding/feature-table-nmds-sample-coordinate.png',
            f'{self.outdir}/feature-embedding/feature-table-nmds-sample-coordinate.pdf',
            f'{self.outdir}/feature-embedding/feature-table-nmds-stress.txt',
        ]:
            self.assertTrue(exists(f))


class TestTSNEProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        TSNEProcess(self.settings).main(
            tsv=f'{self.indir}/feature-table.tsv',
            group_keywords=['H', 'O']
        )
        for f in [
            f'{self.outdir}/feature-embedding/feature-table-tsne-sample-coordinate.tsv',
            f'{self.outdir}/feature-embedding/feature-table-tsne-sample-coordinate.png',
            f'{self.outdir}/feature-embedding/feature-table-tsne-sample-coordinate.pdf',
        ]:
            self.assertTrue(exists(f))
