from .setup import TestCase
from qiime2_pipeline.exporting import ExportFeatureTable, ExportFeatureSequence, \
    ExportTaxonomy, ExportAlignedSequence, ExportTree, ExportBetaDiversity


class TestExportFeatureTable(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ExportFeatureTable(self.settings).main(
            feature_table_qza=f'{self.indir}/feature-table.qza'
        )
        self.assertFileExists(f'{self.workdir}/feature-table.tsv', actual)
        self.assertFileEqual(f'{self.indir}/feature-table.tsv', actual)


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


class TestExportAlignedSequence(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ExportAlignedSequence(self.settings).main(
            aligned_sequence_qza=f'{self.indir}/aligned-sequence.qza'
        )
        expected = f'{self.workdir}/aligned-sequence.fa'
        self.assertFileExists(expected, actual)


class TestExportTree(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ExportTree(self.settings).main(
            tree_qza=f'{self.indir}/tree.qza'
        )
        expected = f'{self.workdir}/tree.nwk'
        self.assertFileExists(expected, actual)


class TestExportBetaDiversity(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ExportBetaDiversity(self.settings).main(
            distance_matrix_qza=f'{self.indir}/distance-matrix.qza'
        )
        expected = f'{self.workdir}/distance-matrix.tsv'
        self.assertFileExists(expected, actual)
