from os.path import exists
from .setup import TestCase
from qiime2_pipeline.lefse import LefSe, OneLefSe


class TestLefSe(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_taxa(self):
        table_tsv_dict = {
            'phylum': f'{self.indir}/taxon-table/phylum-table.tsv',
            'class': f'{self.indir}/taxon-table/class-table.tsv',
            'order': f'{self.indir}/taxon-table/order-table.tsv',
            'family': f'{self.indir}/taxon-table/family-table.tsv',
            'genus': f'{self.indir}/taxon-table/genus-table.tsv',
            'species': f'{self.indir}/taxon-table/species-table.tsv',
        }
        LefSe(self.settings).main(
            table_tsv_dict=table_tsv_dict,
            sample_sheet=f'{self.indir}/taxon-table/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )

    def test_picrust2(self):
        table_tsv_dict = {
            'picrust2-pathway': f'{self.indir}/picrust2/picrust2-pathway-table.tsv',
            'picrust2-EC': f'{self.indir}/picrust2/picrust2-EC-table.tsv',
            'picrust2-KEGG-ortholog': f'{self.indir}/picrust2/picrust2-KEGG-ortholog-table.tsv',
        }
        LefSe(self.settings).main(
            table_tsv_dict=table_tsv_dict,
            sample_sheet=f'{self.indir}/picrust2/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )


class TestOneLefSe(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_genus(self):
        OneLefSe(self.settings).main(
            table_tsv=f'{self.indir}/taxon-table/genus-table.tsv',
            name='genus',
            sample_sheet=f'{self.indir}/taxon-table/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )
        for file in [
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
            f'{self.outdir}/lefse/lefse-picrust2-pathway-features.png',
        ]:
            self.assertTrue(exists(file))
