from .setup import TestCase
from qiime2_pipeline.generate_otu import GenerateOTU, GenerateNanoporeOTU


class TestGenerateOTU(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        clustered_table_qza, clustered_sequence_qza = GenerateOTU(self.settings).main(
            feature_table_qza=f'{self.indir}/dada2-feature-table.qza',
            feature_sequence_qza=f'{self.indir}/dada2-feature-sequence.qza',
            identity=0.975
        )
        for expected, actual in [
            (f'{self.workdir}/otu-feature-table.qza', clustered_table_qza),
            (f'{self.workdir}/otu-feature-sequence.qza', clustered_sequence_qza),
        ]:
            self.assertFileExists(expected, actual)


class TestGenerateNanoporeOTU(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        clustered_table_qza, clustered_sequence_qza = GenerateNanoporeOTU(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet-nanopore.csv',
            fq_dir=f'{self.indir}/fq_dir',
            fq_suffix='_Nanopore.fastq.gz',
            identity=0.975
        )
        for expected, actual in [
            (f'{self.workdir}/otu-feature-table.qza', clustered_table_qza),
            (f'{self.workdir}/otu-feature-sequence.qza', clustered_sequence_qza),
        ]:
            self.assertFileExists(expected, actual)
