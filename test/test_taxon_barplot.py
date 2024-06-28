from os.path import exists
from qiime2_pipeline.taxon_barplot import PlotTaxonBarplots
from .setup import TestCase


class TestPlotTaxonBarplots(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        self.settings.for_publication = True

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
            n_taxa=20,
            sample_sheet=f'{self.indir}/sample-sheet.csv',
        )

        for prefix in ['sample', 'group']:
            for level in [
                'phylum',
                'class',
                'order',
                'family',
                'genus',
                'species',
            ]:
                for ext in ['pdf', 'png', 'tsv']:
                    with self.subTest(prefix=prefix, level=level, ext=ext):
                        self.assertTrue(exists(f'{self.outdir}/taxon-barplot/{prefix}-{level}-barplot.{ext}'))
