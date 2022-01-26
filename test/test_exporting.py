from .setup import TestCase
from qiime2_pipeline.exporting import ExportFeatureTable, ExportFeatureSequence, ExportTaxonomy


class TestExportFeatureTable(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ExportFeatureTable(self.settings).main(
            feature_table_qza=f'{self.indir}/feature-table.qza'
        )
        expected = f'{self.workdir}/feature-table.tsv'
        self.assertFileExists(expected, actual)


class TestExportFeatureSequence(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ExportFeatureSequence(self.settings).main(
            feature_sequence_qza=f'{self.indir}/feature-sequence.qza'
        )
        expected = f'{self.workdir}/feature-sequence.fa'
        self.assertFileExists(expected, actual)


class TestExportTaxonomy(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ExportTaxonomy(self.settings).main(
            taxonomy_qza=f'{self.indir}/taxonomy.qza'
        )
        expected = f'{self.workdir}/taxonomy.tsv'
        self.assertFileExists(expected, actual)
