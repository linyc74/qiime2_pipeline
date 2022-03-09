import pandas as pd
from qiime2_pipeline.heatmap import Clustermap, FilterByCumulativeReads, PlotHeatmaps
from .setup import TestCase


class TestPlotHeatmaps(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        tsvs = [
            f'{self.indir}/class.tsv',
            f'{self.indir}/species.tsv',
        ]
        PlotHeatmaps(self.settings).main(
            tsvs=tsvs,
            heatmap_read_fraction=0.9
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


class TestClustermap(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        Clustermap(self.settings).main(
            data=pd.read_csv(f'{self.indir}/filtered.tsv', sep='\t', index_col=0),
            output_prefix=f'{self.outdir}/clustermap'
        )
