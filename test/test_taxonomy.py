from .setup import TestCase
from qiime2_pipeline.taxonomy import Taxonomy, MergeForwardReverseTaxonomy


class TestTaxonomy(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_nb(self):
        actual = Taxonomy(self.settings).main(
            representative_seq_qza=f'{self.indir}/dada2-feature-sequence.qza',
            feature_classifier='nb',
            nb_classifier_qza=f'{self.indir}/gg-13-8-99-515-806-nb-classifier.qza',
            classifier_reads_per_batch=0,
            reference_sequence_qza=None,
            reference_taxonomy_qza=None
        )
        expected = f'{self.workdir}/taxonomy-merged.qza'
        self.assertFileExists(expected, actual)

    def test_vsearch(self):
        actual = Taxonomy(self.settings).main(
            representative_seq_qza=f'{self.indir}/dada2-feature-sequence.qza',
            feature_classifier='vsearch',
            nb_classifier_qza=None,
            classifier_reads_per_batch=0,
            reference_sequence_qza=f'{self.indir}/silva-138-99-seqs.qza',
            reference_taxonomy_qza=f'{self.indir}/silva-138-99-tax.qza'
        )
        expected = f'{self.workdir}/taxonomy-vsearch.qza'
        self.assertFileExists(expected, actual)


class TestMergeForwardReverseTaxonomy(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = MergeForwardReverseTaxonomy(self.settings).main(
            forward_taxonomy_qza=f'{self.indir}/taxonomy-same.qza',
            reverse_taxonomy_qza=f'{self.indir}/taxonomy-reverse-complement.qza',
        )
        expected = f'{self.workdir}/taxonomy-merged.qza'
        self.assertFileExists(expected, actual)
