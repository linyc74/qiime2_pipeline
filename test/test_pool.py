from .setup import TestCase
from qiime2_pipeline.pool import Pool, BatchPool


class TestPool(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = Pool(self.settings).main(
            fq1=f'{self.indir}/R1.fastq.gz',
            fq2=f'{self.indir}/R2.fastq.gz'
        )
        expected = f'{self.workdir}/pool.fq'
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
        self.assertEqual('.fq', fq_suffix)
