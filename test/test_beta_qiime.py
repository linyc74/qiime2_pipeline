from .setup import TestCase
from qiime2_pipeline.beta_qiime import QiimeBetaDiversity, RunAllBetaMetricsToTsvs


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


class TestRunAllBetaMetricsToTsvs(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = RunAllBetaMetricsToTsvs(self.settings).main(
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
