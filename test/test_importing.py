from .setup import TestCase
from qiime2_pipeline.exporting import ExportFeatureTable, ExportFeatureSequence
from qiime2_pipeline.importing import ImportFeatureTable, ImportFeatureSequence, ImportTaxonomy, \
    ImportPairedEndFastq, ImportSingleEndFastq


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


class TestImportTaxonomy(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ImportTaxonomy(self.settings).main(
            taxonomy_tsv=f'{self.indir}/taxonomy.tsv')
        expected = f'{self.workdir}/taxonomy.qza'
        self.assertFileExists(expected, actual)


class TestImportPairedEndFastq(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ImportPairedEndFastq(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz')
        expected = f'{self.workdir}/paired-end-demultiplexed.qza'
        self.assertFileExists(expected, actual)


class TestImportSingleEndFastq(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ImportSingleEndFastq(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            fq_dir=f'{self.indir}/concat_fastqs',
            fq_suffix='.fq'
        )
        expected = f'{self.workdir}/single-end-demultiplexed.qza'
        self.assertFileExists(expected, actual)
