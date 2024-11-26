from .setup import TestCase
from qiime2_pipeline.generate_otu import GenerateOTUFromASV


class TestGenerateOTUFromASV(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        clustered_table_qza, clustered_sequence_qza = GenerateOTUFromASV(self.settings).main(
            feature_table_qza=f'{self.indir}/dada2-feature-table.qza',
            feature_sequence_qza=f'{self.indir}/dada2-feature-sequence.qza',
            identity=0.975
        )
        for expected, actual in [
            (f'{self.workdir}/vsearch-feature-table.qza', clustered_table_qza),
            (f'{self.workdir}/vsearch-feature-sequence.qza', clustered_sequence_qza),
        ]:
            self.assertFileExists(expected, actual)



