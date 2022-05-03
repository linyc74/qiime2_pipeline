from qiime2_pipeline.phylogeny import Phylogeny, DrawTree
from .setup import TestCase


class TestPhylogeny(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test_main(self):
        actual = Phylogeny(self.settings).main(
            seq_qza=f'{self.indir}/labeled-feature-sequence.qza')
        expected = f'{self.workdir}/fasttree-rooted.qza'
        self.assertFileExists(expected, actual)


class TestDrawTree(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test_main(self):
        actual = DrawTree(self.settings).main(
            nwk=f'{self.indir}/fasttree-rooted.nwk'
        )
