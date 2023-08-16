from qiime2_pipeline.differential_abundance import DifferentialAbundance, OneTaxonLevelDifferentialAbundance
from .setup import TestCase


class TestDifferentialAbundance(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        taxon_table_tsv_dict = {
            'phylum': f'{self.indir}/phylum-table.tsv',
            'class': f'{self.indir}/class-table.tsv',
            'genus': f'{self.indir}/genus-table.tsv',
        }
        DifferentialAbundance(self.settings).main(
            taxon_table_tsv_dict=taxon_table_tsv_dict,
            group_keywords=['H', 'O']
        )


class TestOneTaxonLevelDifferentialAbundance(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        OneTaxonLevelDifferentialAbundance(self.settings).main(
            taxon_level='genus',
            taxon_tsv=f'{self.indir}/genus-table.tsv',
            group_keywords=['H', 'O', 'A']
        )

    def test_error_separator(self):
        OneTaxonLevelDifferentialAbundance(self.settings).main(
            taxon_level='species',
            taxon_tsv=f'{self.indir}/error-separator-table.tsv',
            group_keywords=['CA', 'N']
        )
