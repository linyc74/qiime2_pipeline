from qiime2_pipeline.phylogeny import Phylogeny, DrawTree, leaf_name_to_bacterial_phylum
from .setup import TestCase


class TestPhylogeny(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = Phylogeny(self.settings).main(
            seq_qza=f'{self.indir}/labeled-feature-sequence.qza')
        expected = f'{self.workdir}/fasttree-rooted.qza'
        self.assertFileExists(expected, actual)


class TestDrawTree(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        DrawTree(self.settings).main(
            nwk=f'{self.indir}/fasttree-rooted-big.nwk',
            dstdir=self.outdir)

    def test_leaf_name_to_bacterial_phylum(self):
        func = leaf_name_to_bacterial_phylum

        leaf_name = 'ASV_0053; d__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales'
        self.assertEqual('Bacteroidetes', func(leaf_name))

        leaf_name = 'ASV_0053; d__Bacteria; p__Bacteroidetes'
        self.assertEqual('Bacteroidetes', func(leaf_name))

        leaf_name = 'ASV_0053; p__Bacteroidetes'
        self.assertEqual('Unassigned', func(leaf_name))

        leaf_name = 'ASV_0053; d__Bacteria'
        self.assertEqual('Unassigned', func(leaf_name))
