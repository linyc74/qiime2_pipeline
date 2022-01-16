from .setup import TestCase
from qiime2_pipeline.denoise import Dada2Paired
from qiime2_pipeline.trimming import TrimGalore, BatchTrimGalore
from qiime2_pipeline.importing import ImportPairedEndFastq
from qiime2_pipeline.concat import Concat, BatchConcat


class MyTest(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def __test_trim_galore(self):
        trimmed_fq1, trimmed_fq2 = TrimGalore(self.settings).main(
            fq1=f'{self.indir}/R1.fastq.gz',
            fq2=f'{self.indir}/R2.fastq.gz',
        )
        for expected, actual in [
            (f'{self.workdir}/R1_val_1.fq.gz', trimmed_fq1),
            (f'{self.workdir}/R2_val_2.fq.gz', trimmed_fq2),
        ]:
            self.assertFileExists(expected, actual)

    def __test_batch_trim_galore(self):
        actual = BatchTrimGalore(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz')
        expected = f'{self.workdir}/trimmed_fastqs'
        self.assertFileExists(expected, actual)

    def __test_concat(self):
        actual = Concat(self.settings).main(
            fq1=f'{self.indir}/R1.fastq.gz',
            fq2=f'{self.indir}/R2.fastq.gz'
        )
        expected = f'{self.workdir}/concat.fq'
        self.assertFileExists(expected, actual)

    def __test_batch_concat(self):
        fq_dir, fq_suffix = BatchConcat(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz')
        self.assertFileExists(f'{self.workdir}/concat_fastqs', fq_dir)
        self.assertEqual('.fq', fq_suffix)

    def __test_import_paired_end_fastq(self):
        actual = ImportPairedEndFastq(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz')
        expected = f'{self.workdir}/paired-end-demultiplexed.qza'
        self.assertFileExists(expected, actual)

    def __test_dada2_paired(self):
        Dada2Paired(self.settings).main(
            demultiplexed_seq_qza=f'{self.indir}/paired-end-demultiplexed.qza'
        )
