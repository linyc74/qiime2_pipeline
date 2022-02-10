from .setup import TestCase
from qiime2_pipeline.generate_asv import GenerateASVConcatPairedEnd, GenerateASVPoolPairedEnd, GenerateASVSingleEnd


class TestGenerateASVConcatPairedEnd(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        feature_table_qza, feature_sequence_qza = GenerateASVConcatPairedEnd(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz',
            clip_r1_5_prime=17,
            clip_r2_5_prime=0
        )
        for expected, actual in [
            (f'{self.workdir}/dada2-feature-table.qza', feature_table_qza),
            (f'{self.workdir}/dada2-feature-sequence.qza', feature_sequence_qza),
        ]:
            self.assertFileExists(expected, actual)


class TestGenerateASVPoolPairedEnd(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        feature_table_qza, feature_sequence_qza = GenerateASVPoolPairedEnd(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz',
            clip_r1_5_prime=17,
            clip_r2_5_prime=0
        )
        for expected, actual in [
            (f'{self.workdir}/dada2-feature-table.qza', feature_table_qza),
            (f'{self.workdir}/dada2-feature-sequence.qza', feature_sequence_qza),
        ]:
            self.assertFileExists(expected, actual)


class TestGenerateASVSingleEnd(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        feature_table_qza, feature_sequence_qza = GenerateASVSingleEnd(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq_suffix='_L001_R1_001.fastq.gz',
            clip_5_prime=17
        )
        for expected, actual in [
            (f'{self.workdir}/dada2-feature-table.qza', feature_table_qza),
            (f'{self.workdir}/dada2-feature-sequence.qza', feature_sequence_qza),
        ]:
            self.assertFileExists(expected, actual)
