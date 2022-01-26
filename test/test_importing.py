from .setup import TestCase
from qiime2_pipeline.importing import ImportFeatureTable, ImportFeatureSequence
from qiime2_pipeline.exporting import ExportFeatureTable, ExportFeatureSequence


class TestImportFeatureTable(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ImportFeatureTable(self.settings).main(
            feature_table_tsv=f'{self.indir}/feature-table.tsv')
        expected = f'{self.workdir}/feature-table.qza'
        self.assertFileExists(expected, actual)

    def test_data_integrity(self):
        qza = ImportFeatureTable(self.settings).main(
            feature_table_tsv=f'{self.indir}/feature-table.tsv')
        actual = ExportFeatureTable(self.settings).main(
            feature_table_qza=qza)
        expected = f'{self.indir}/feature-table-from-qza.tsv'
        self.assertFileEqual(expected, actual)


class TestImportFeatureSequence(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ImportFeatureSequence(self.settings).main(
            feature_sequence_fa=f'{self.indir}/feature-sequence.fa')
        expected = f'{self.workdir}/feature-sequence.qza'
        self.assertFileExists(expected, actual)

    def test_data_integrity(self):
        input_fa = f'{self.indir}/feature-sequence.fa'
        qza = ImportFeatureSequence(self.settings).main(feature_sequence_fa=input_fa)
        output_fa = ExportFeatureSequence(self.settings).main(feature_sequence_qza=qza)
        self.assertFileEqual(input_fa, output_fa)
