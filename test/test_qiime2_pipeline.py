from qiime2_pipeline.qiime2_pipeline import Qiime2Pipeline
from .setup import TestCase


class TestQiime2Pipeline(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_paired_end(self):
        Qiime2Pipeline(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz',
            nb_classifier_qza=f'{self.indir}/gg-13-8-99-515-806-nb-classifier.qza',
            paired_end_mode='pool',
            group_keywords=['H', 'O'],
            otu_identity=0.97,
            skip_otu=False,
            classifier_reads_per_batch=0,
            alpha_metrics=['shannon', 'observed_features']
        )

    def test_single_end(self):
        Qiime2Pipeline(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix=None,
            nb_classifier_qza=f'{self.indir}/gg-13-8-99-515-806-nb-classifier.qza',
            paired_end_mode='pool',
            group_keywords=['H', 'O'],
            otu_identity=0.97,
            skip_otu=False,
            classifier_reads_per_batch=0,
            alpha_metrics=['shannon', 'observed_features']
        )
