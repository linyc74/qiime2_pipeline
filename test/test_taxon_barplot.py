from os.path import exists
from qiime2_pipeline.taxon_barplot import PlotTaxonBarplots
from .setup import TestCase


class TestPlotTaxonBarplots(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        taxon_table_tsv_dict = {
            'phylum': f'{self.indir}/phylum-table.tsv',
            'class': f'{self.indir}/class-table.tsv',
            'order': f'{self.indir}/order-table.tsv',
            'family': f'{self.indir}/family-table.tsv',
            'genus': f'{self.indir}/genus-table.tsv',
            'species': f'{self.indir}/species-table.tsv',
        }

        PlotTaxonBarplots(self.settings).main(
            taxon_table_tsv_dict=taxon_table_tsv_dict,
            n_taxa=20
        )

        for level in [
            'phylum',
            'class',
            'order',
            'family',
            'genus',
            'species',
        ]:
            for ext in ['png', 'tsv']:
                fpath = f'{self.outdir}/taxon-barplot/{level}-barplot.{ext}'
                self.assertTrue(exists(fpath))
