from os.path import exists
from .setup import TestCase
from qiime2_pipeline.lefse import LefSe, LefSeOneTaxonLevel, InsertGroupRow


class TestLefSe(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        taxon_table_tsv_dict = {
            'phylum': f'{self.indir}/taxon-table/phylum-table.tsv',
            'class': f'{self.indir}/taxon-table/class-table.tsv',
            'order': f'{self.indir}/taxon-table/order-table.tsv',
            'family': f'{self.indir}/taxon-table/family-table.tsv',
            'genus': f'{self.indir}/taxon-table/genus-table.tsv',
            'species': f'{self.indir}/taxon-table/species-table.tsv',
        }
        LefSe(self.settings).main(
            taxon_table_tsv_dict=taxon_table_tsv_dict,
            group_keywords=['H', 'O']
        )


class TestLefSeOneTaxonLevel(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        LefSeOneTaxonLevel(self.settings).main(
            taxon_table_tsv=f'{self.indir}/genus-table.tsv',
            taxon_level='genus',
            group_keywords=['H', 'O']
        )
        for file in [
            f'{self.outdir}/lefse/lefse-genus-result.txt',
            f'{self.outdir}/lefse/lefse-genus-features.png',
            f'{self.outdir}/lefse/lefse-genus-cladogram.png',
        ]:
            self.assertTrue(exists(file))


class TestInsertGroupRow(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = InsertGroupRow(self.settings).main(
            taxon_table_tsv=f'{self.indir}/genus-table.tsv',
            group_keywords=['H']
        )
        expected = f'{self.indir}/genus-table-grouped.tsv'
        self.assertFileEqual(expected, actual)
