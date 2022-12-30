from .setup import TestCase
from qiime2_pipeline.raw_read_counts import RawReadCounts


class TestReadCounts(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_single_end(self):
        RawReadCounts(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix=None,
        )

    def test_paired_end(self):
        RawReadCounts(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz',
        )
