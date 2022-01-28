from .setup import TestCase
from os.path import exists
from qiime2_pipeline.dim_reduction import PCoA, NMDS, TSNE


class TestPCoA(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        PCoA(self.settings).main(
            distance_matrix_tsv=f'{self.indir}/distance-matrix.tsv'
        )
        for f in [
            f'{self.outdir}/PCoA/distance-matrix-pcoa-sample-coordinate.tsv',
            f'{self.outdir}/PCoA/distance-matrix-pcoa-proportion-explained.tsv',
        ]:
            self.assertTrue(exists(f))


class TestNMDS(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        NMDS(self.settings).main(
            distance_matrix_tsv=f'{self.indir}/distance-matrix.tsv'
        )
        for f in [
            f'{self.outdir}/NMDS/distance-matrix-nmds-sample-coordinate.tsv',
            f'{self.outdir}/NMDS/distance-matrix-nmds-stress.txt',
        ]:
            self.assertTrue(exists(f))


class TestTSNE(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        TSNE(self.settings).main(
            distance_matrix_tsv=f'{self.indir}/distance-matrix.tsv'
        )
        for f in [
            f'{self.outdir}/t-SNE/distance-matrix-tsne-sample-coordinate.tsv',
        ]:
            self.assertTrue(exists(f))
