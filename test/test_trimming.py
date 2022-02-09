from .setup import TestCase
from qiime2_pipeline.trimming import TrimGalorePairedEnd, BatchTrimGalorePairedEnd, \
    TrimGaloreSingleEnd, BatchTrimGaloreSingleEnd


class TestTrimGalorePairedEnd(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        trimmed_fq1, trimmed_fq2 = TrimGalorePairedEnd(self.settings).main(
            fq1=f'{self.indir}/R1.fastq.gz',
            fq2=f'{self.indir}/R2.fastq.gz',
        )
        for expected, actual in [
            (f'{self.workdir}/R1_val_1.fq.gz', trimmed_fq1),
            (f'{self.workdir}/R2_val_2.fq.gz', trimmed_fq2),
        ]:
            self.assertFileExists(expected, actual)


class TestBatchTrimGalorePairedEnd(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = BatchTrimGalorePairedEnd(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz')
        expected = f'{self.workdir}/trimmed_fastqs'
        self.assertFileExists(expected, actual)


class TestTrimGaloreSingleEnd(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        trimmed_fq = TrimGaloreSingleEnd(self.settings).main(
            fq=f'{self.indir}/R1.fastq.gz',
        )
        self.assertFileExists(f'{self.workdir}/R1_trimmed.fq.gz', trimmed_fq)


class TestBatchTrimGaloreSingleEnd(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = BatchTrimGaloreSingleEnd(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq_suffix='_L001_R2_001.fastq.gz')
        expected = f'{self.workdir}/trimmed_fastqs'
        self.assertFileExists(expected, actual)
