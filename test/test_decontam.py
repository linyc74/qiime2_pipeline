from .setup import TestCase
from qiime2_pipeline.decontam import Decontam


class TestDecontam(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        decontam_table_qza, decontam_sequence_qza = Decontam(self.settings).main(
            feature_table_qza=f'{self.indir}/feature-table.qza',
            feature_sequence_qza=f'{self.indir}/feature-sequence.qza',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            dna_concentration_column='DNA conc. (ng/µL)',
            decontam_threshold=0.1,
        )
        for expected, actual in [
            (f'{self.workdir}/decontam-feature-table.qza', decontam_table_qza),
            (f'{self.workdir}/decontam-feature-sequence.qza', decontam_sequence_qza),
        ]:
            self.assertFileExists(expected, actual)
