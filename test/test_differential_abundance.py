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
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0), (0.1, 0.9, 0.5, 1.0), (0.6, 0.2, 0.4, 1.0)],
            p_value=0.05
        )

    def test_number_groups(self):
        taxon_table_tsv_dict = {
            'phylum': f'{self.indir}/phylum-table.tsv',
            'class': f'{self.indir}/class-table.tsv',
            'genus': f'{self.indir}/genus-table.tsv',
        }
        DifferentialAbundance(self.settings).main(
            taxon_table_tsv_dict=taxon_table_tsv_dict,
            sample_sheet=f'{self.indir}/sample-sheet-number-groups.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0), (0.1, 0.9, 0.5, 1.0), (0.6, 0.2, 0.4, 1.0)],
            p_value=0.05
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
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0), (0.1, 0.9, 0.5, 1.0), (0.6, 0.2, 0.4, 1.0)],
            p_value=0.05
        )

    def test_error_separator(self):
        OneTaxonLevelDifferentialAbundance(self.settings).main(
            taxon_level='genus',
            taxon_tsv=f'{self.indir}/error-separator-table.tsv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)],
            p_value=0.05
        )
