from os.path import exists
from .setup import TestCase
from qiime2_pipeline.lefse import LefSe, OneLefSe, InsertGroupRow, AddTaxonLevelPrefixes


class TestLefSe(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        indir = f'{self.indir}/mock-taxon-table'
        taxon_table_tsv_dict = {
            'phylum': f'{indir}/phylum-table.tsv',
            'class': f'{indir}/class-table.tsv',
            'order': f'{indir}/order-table.tsv',
            'family': f'{indir}/family-table.tsv',
            'genus': f'{indir}/genus-table.tsv',
            'species': f'{indir}/species-table.tsv',
        }
        LefSe(self.settings).main(
            table_tsv_dict=taxon_table_tsv_dict,
            sample_sheet=f'{indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )


class TestOneLefSe(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_genus(self):
        OneLefSe(self.settings).main(
            table_tsv=f'{self.indir}/real-taxon-table/genus-table.tsv',
            name='genus',
            sample_sheet=f'{self.indir}/real-taxon-table/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )
        for file in [
            f'{self.outdir}/lefse/lefse-genus-result.txt',
            f'{self.outdir}/lefse/lefse-genus-features.png',
            f'{self.outdir}/lefse/lefse-genus-cladogram.png',
        ]:
            self.assertTrue(exists(file))

    def test_picrust2(self):
        OneLefSe(self.settings).main(
            table_tsv=f'{self.indir}/picrust2/picrust2-pathway-table.tsv',
            name='picrust2-pathway',
            sample_sheet=f'{self.indir}/picrust2/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )
        for file in [
            f'{self.outdir}/lefse/lefse-picrust2-pathway-result.txt',
            f'{self.outdir}/lefse/lefse-picrust2-pathway-features.png',
        ]:
            self.assertTrue(exists(file))


class TestAddTaxonLevelPrefixes(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = AddTaxonLevelPrefixes(self.settings).main(
            taxon_table_tsv=f'{self.indir}/mock-taxon-table/genus-table.tsv',
        )
        expected = f'{self.indir}/mock-taxon-table/genus-table-relabeled.tsv'
        self.assertFileEqual(expected, actual)

    def test_add_level_prefixes(self):
        actual = AddTaxonLevelPrefixes(self.settings).add_level_prefixes(
            s='A|B|C|D|E|F|G'
        )
        expected = 'domain__A|phylum__B|class__C|order__D|family__E|genus__F|species__G'
        self.assertEqual(actual, expected)


class TestInsertGroupRow(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = InsertGroupRow(self.settings).main(
            table_tsv=f'{self.indir}/mock-taxon-table/genus-table.tsv',
            sample_sheet=f'{self.indir}/mock-taxon-table/sample-sheet.csv',
        )
        expected = f'{self.indir}/mock-taxon-table/genus-table-grouped.tsv'
        self.assertFileEqual(expected, actual)
