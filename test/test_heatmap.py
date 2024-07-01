import pandas as pd
from qiime2_pipeline.heatmap import PlotHeatmaps, PlotOneHeatmap, FilterByCumulativeReads
from .setup import TestCase


class TestPlotHeatmaps(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        tsvs = [
            f'{self.indir}/phylum-table.tsv',
            f'{self.indir}/class-table.tsv',
            f'{self.indir}/order-table.tsv',
            f'{self.indir}/family-table.tsv',
            f'{self.indir}/genus-table.tsv',
            f'{self.indir}/species-table.tsv',
        ]
        PlotHeatmaps(self.settings).main(
            tsvs=tsvs,
            heatmap_read_fraction=0.9,
            sample_sheet=f'{self.indir}/sample-sheet.csv',
        )


class TestPlotOneHeatmap(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        self.settings.for_publication = False
        PlotOneHeatmap(self.settings).main(
            tsv=f'{self.indir}/346-samples-species-table.tsv',
            heatmap_read_fraction=0.95,
            sample_sheet=f'{self.indir}/346-samples-sample-sheet.csv',
            dstdir=self.outdir
        )


class TestFilterByCumulativeReads(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = FilterByCumulativeReads(self.settings).main(
            df=pd.read_csv(f'{self.indir}/unfiltered.tsv', sep='\t', index_col=0),
            heatmap_read_fraction=0.8
        )
        expected = pd.read_csv(f'{self.indir}/filtered.tsv', sep='\t', index_col=0)
        self.assertDataFrameEqual(expected, actual)
