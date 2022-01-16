from .setup import TestCase
from qiime2_pipeline.denoise import Dada2Paired
from qiime2_pipeline.trimming import CutadaptTrimPaired
from qiime2_pipeline.importing import ImportPairedEndFastq


class MyTest(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def __test_import_and_trimming(self):
        untrimmed_reads_qza = ImportPairedEndFastq(self.settings).main(
            fq_dir=f'{self.indir}/data',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz')

        trimmed_reads_qza = CutadaptTrimPaired(self.settings).main(
            untrimmed_reads_qza=untrimmed_reads_qza)

    def __test_dada2_paired(self):
        Dada2Paired(self.settings).main(
            demultiplexed_seq_qza=f'{self.indir}/paired-end-demultiplexed.qza'
        )
