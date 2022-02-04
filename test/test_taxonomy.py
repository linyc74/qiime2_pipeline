from .setup import TestCase
from qiime2_pipeline.taxonomy import Taxonomy, Classify, MergeForwardReverseTaxonomy


class TestTaxonomy(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = Taxonomy(self.settings).main(
            representative_seq_qza=f'{self.indir}/dada2-feature-sequence.qza',
            nb_classifier_qza=f'{self.indir}/gg-13-8-99-515-806-nb-classifier.qza'
        )
        expected = f'{self.workdir}/taxonomy-merged.qza'
        self.assertFileExists(expected, actual)


class TestClassify(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def __test_main(self):
        actual = Classify(self.settings).main(
            representative_seq_qza=f'{self.indir}/dada2-feature-sequence.qza',
            nb_classifier_qza=f'{self.indir}/gg-13-8-99-515-806-nb-classifier.qza',
            read_orientation='same'
        )
        expected = f'{self.workdir}/taxonomy-same.qza'
        self.assertFileExists(expected, actual)


class TestMergeForwardReverseTaxonomy(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def __test_main(self):
        actual = MergeForwardReverseTaxonomy(self.settings).main(
            forward_taxonomy_qza=f'{self.indir}/taxonomy-same.qza',
            reverse_taxonomy_qza=f'{self.indir}/taxonomy-reverse-complement.qza',
        )
        expected = f'{self.workdir}/taxonomy-merged.qza'
        self.assertFileExists(expected, actual)
