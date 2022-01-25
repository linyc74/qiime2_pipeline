from .setup import TestCase
from qiime2_pipeline.labeling import FeatureLabeling, ImportFeatureTable, ImportFeatureSequence


class TestFeatureLabeling(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        table_qza, sequence_qza = FeatureLabeling(self.settings).main(
            taxonomy_qza=f'{self.indir}/taxonomy.qza',
            feature_table_qza=f'{self.indir}/dada2-feature-table.qza',
            feature_sequence_qza=f'{self.indir}/dada2-feature-sequence.qza'
        )
        for expected, actual in [
            (f'{self.workdir}/labeled-feature-sequence.qza', sequence_qza),
            (f'{self.workdir}/labeled-feature-table.qza', table_qza)
        ]:
            self.assertFileExists(expected, actual)


class TestImportFeatureTable(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ImportFeatureTable(self.settings).main(
            feature_table_tsv=f'{self.indir}/labeled-feature-table.tsv')
        expected = f'{self.workdir}/labeled-feature-table.qza'
        self.assertFileExists(expected, actual)


class TestImportFeatureSequence(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ImportFeatureSequence(self.settings).main(
            feature_sequence_fa=f'{self.indir}/labeled-feature-sequence.fa')
        expected = f'{self.workdir}/labeled-feature-sequence.qza'
        self.assertFileExists(expected, actual)
