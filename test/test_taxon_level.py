import pandas as pd
from os.path import exists
from qiime2_pipeline.taxon_level import TaxonLevel, CollapseTaxon, feature_label_to_taxon, trim_taxon
from .setup import TestCase


class TestTaxonLevel(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = TaxonLevel(self.settings).main(
            labeled_feature_table_tsv=f'{self.indir}/labeled-feature-table.tsv'
        )
        expected = {
            'phylum': f'{self.outdir}/taxon-level/phylum.tsv',
            'class': f'{self.outdir}/taxon-level/class.tsv',
            'order': f'{self.outdir}/taxon-level/order.tsv',
            'family': f'{self.outdir}/taxon-level/family.tsv',
            'genus': f'{self.outdir}/taxon-level/genus.tsv',
            'species': f'{self.outdir}/taxon-level/species.tsv',
        }
        self.assertDictEqual(expected, actual)
        for path in expected.values():
            self.assertTrue(exists(path))


class TestCollapseTaxon(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = CollapseTaxon(self.settings).main(
            df=pd.read_csv(f'{self.indir}/species.tsv', sep='\t', index_col=0),
            level='class'
        )
        expected = pd.read_csv(f'{self.indir}/class.tsv', sep='\t', index_col=0)
        self.assertDataFrameEqual(expected, actual)


class TestFunction(TestCase):

    def test_feature_label_to_taxon(self):
        actual = feature_label_to_taxon(
            s='X; x__AAA; x__BBB; x__CCC; x__DDD; x__EEE; x__FFF; x__GGG'
        )
        expected = 'AAA|BBB|CCC|DDD|EEE|FFF|GGG'
        self.assertEqual(expected, actual)

    def test_trim_taxon(self):
        actual = trim_taxon(
            taxon='1|2|3|4|5|6|7',
            int_level=5
        )
        expected = '1|2|3|4|5'
        self.assertEqual(expected, actual)
