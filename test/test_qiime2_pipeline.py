from qiime2_pipeline.qiime2_pipeline import Qiime2Pipeline
from .setup import TestCase


class TestQiime2Pipeline(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_paired_end(self):
        Qiime2Pipeline(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz',
            pacbio=False,
            paired_end_mode='pool',
            otu_identity=0.97,
            skip_otu=False,
            feature_classifier='nb',
            nb_classifier_qza=f'{self.indir}/gg-13-8-99-515-806-nb-classifier.qza',
            classifier_reads_per_batch=0,
            reference_sequence_qza=None,
            reference_taxonomy_qza=None,
            alpha_metrics=['shannon', 'observed_features'],
            clip_r1_5_prime=17,
            clip_r2_5_prime=0,
            max_expected_error_bases=2.0,
            heatmap_read_fraction=0.99,
            n_taxa_barplot=20,
            beta_diversity_feature_level='feature',
            colormap='gold,limegreen',
            invert_colors=True,
            skip_differential_abundance=False,
            run_picrust2=False,
        )

    def test_single_end(self):
        Qiime2Pipeline(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix=None,
            pacbio=False,
            paired_end_mode='pool',
            otu_identity=0.97,
            skip_otu=False,
            feature_classifier='nb',
            nb_classifier_qza=f'{self.indir}/gg-13-8-99-515-806-nb-classifier.qza',
            classifier_reads_per_batch=0,
            reference_sequence_qza=None,
            reference_taxonomy_qza=None,
            alpha_metrics=['shannon', 'observed_features'],
            clip_r1_5_prime=17,
            clip_r2_5_prime=0,
            max_expected_error_bases=2.0,
            heatmap_read_fraction=0.99,
            n_taxa_barplot=20,
            beta_diversity_feature_level='feature',
            colormap='Set1',
            invert_colors=True,
            skip_differential_abundance=False,
            run_picrust2=False,
        )

    def test_fungi(self):
        Qiime2Pipeline(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz',
            pacbio=False,
            paired_end_mode='pool',
            otu_identity=0.97,
            skip_otu=True,
            feature_classifier='nb',
            nb_classifier_qza=f'{self.indir}/unite_ver8_dynamic_10.05.2021-Q2-2022.8.qza',
            classifier_reads_per_batch=0,
            reference_sequence_qza=None,
            reference_taxonomy_qza=None,
            alpha_metrics=['shannon', 'observed_features'],
            clip_r1_5_prime=17,
            clip_r2_5_prime=0,
            max_expected_error_bases=2.0,
            heatmap_read_fraction=0.99,
            n_taxa_barplot=20,
            beta_diversity_feature_level='species',
            colormap='gold,limegreen',
            invert_colors=True,
            skip_differential_abundance=False,
            run_picrust2=False,
        )

    def test_pacbio(self):
        Qiime2Pipeline(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet-pacbio.csv',
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='.fastq.gz',
            fq2_suffix=None,
            pacbio=True,
            paired_end_mode='merge',
            otu_identity=0.97,
            skip_otu=False,
            feature_classifier='nb',
            nb_classifier_qza=f'{self.indir}/gg-13-8-99-nb-classifier.qza',
            classifier_reads_per_batch=0,
            reference_sequence_qza=None,
            reference_taxonomy_qza=None,
            alpha_metrics=['shannon', 'observed_features'],
            clip_r1_5_prime=0,
            clip_r2_5_prime=0,
            max_expected_error_bases=8.0,
            heatmap_read_fraction=0.99,
            n_taxa_barplot=20,
            beta_diversity_feature_level='feature',
            colormap='Set1',
            invert_colors=True,
            skip_differential_abundance=False,
            run_picrust2=False,
        )