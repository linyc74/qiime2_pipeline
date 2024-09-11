from .setup import TestCase
from qiime2_pipeline.generate_asv import GenerateASV, BatchPool


class TestGenerateASV(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_pacbio(self):
        feature_table_qza, feature_sequence_qza = GenerateASV(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet-pacbio.csv',
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='.fastq.gz',
            fq2_suffix=None,
            pacbio=True,
            paired_end_mode='merge',  # not used anyway
            clip_r1_5_prime=0,
            clip_r2_5_prime=0,
            max_expected_error_bases=8.,
        )
        for expected, actual in [
            (f'{self.workdir}/dada2-feature-table.qza', feature_table_qza),
            (f'{self.workdir}/dada2-feature-sequence.qza', feature_sequence_qza),
        ]:
            self.assertFileExists(expected, actual)

    def test_single_end(self):
        feature_table_qza, feature_sequence_qza = GenerateASV(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix=None,
            pacbio=False,
            paired_end_mode='merge',  # not used anyway
            clip_r1_5_prime=17,
            clip_r2_5_prime=0,
            max_expected_error_bases=2.,
        )
        for expected, actual in [
            (f'{self.workdir}/dada2-feature-table.qza', feature_table_qza),
            (f'{self.workdir}/dada2-feature-sequence.qza', feature_sequence_qza),
        ]:
            self.assertFileExists(expected, actual)

    def test_paired_end_pool(self):
        feature_table_qza, feature_sequence_qza = GenerateASV(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz',
            pacbio=False,
            paired_end_mode='pool',
            clip_r1_5_prime=17,
            clip_r2_5_prime=0,
            max_expected_error_bases=2.,
        )
        for expected, actual in [
            (f'{self.workdir}/dada2-feature-table.qza', feature_table_qza),
            (f'{self.workdir}/dada2-feature-sequence.qza', feature_sequence_qza),
        ]:
            self.assertFileExists(expected, actual)

    def test_paired_end_merge(self):
        feature_table_qza, feature_sequence_qza = GenerateASV(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz',
            pacbio=False,
            paired_end_mode='merge',
            clip_r1_5_prime=17,
            clip_r2_5_prime=0,
            max_expected_error_bases=2.,
        )
        for expected, actual in [
            (f'{self.workdir}/dada2-feature-table.qza', feature_table_qza),
            (f'{self.workdir}/dada2-feature-sequence.qza', feature_sequence_qza),
        ]:
            self.assertFileExists(expected, actual)


class TestBatchPool(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        fq_dir, fq_suffix = BatchPool(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz')
        self.assertFileExists(f'{self.workdir}/pool_fastqs', fq_dir)
        self.assertEqual('.fastq.gz', fq_suffix)
