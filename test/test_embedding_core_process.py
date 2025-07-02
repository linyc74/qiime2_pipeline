import pandas as pd
from qiime2_pipeline.embedding_core_process import PCACore, ScatterPlot
from .setup import TestCase


class TestCores(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)
        self.feature_table_df = read_tsv(f'{self.indir}/feature-table.tsv')

    def tearDown(self):
        self.tear_down()

    def test_pca_core(self):
        sample_coordinate_df, proportion_explained_series = PCACore(self.settings).main(
            df=self.feature_table_df,
            data_structure='row_features'
        )
        expected = read_tsv(f'{self.indir}/pca_sample_coordinate_df.tsv')
        self.assertDataFrameEqual(expected, sample_coordinate_df)


class TestScatterPlot(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_short_name(self):
        ScatterPlot(self.settings).main(
            sample_coordinate_df=read_tsv(f'{self.indir}/sample-coordinate-short-name.tsv'),
            x_column='PC1',
            y_column='PC2',
            hue_column='Group',
            colors=['#1f77b4', '#ff7f0e', '#2ca02c'],
            output_prefix=f'{self.outdir}/scatterplot')

    def test_long_name(self):
        ScatterPlot(self.settings).main(
            sample_coordinate_df=read_tsv(f'{self.indir}/sample-coordinate-long-name.tsv'),
            x_column='PC1',
            y_column='PC2',
            hue_column='Group',
            colors=['#1f77b4', '#ff7f0e', '#2ca02c'],
            output_prefix=f'{self.outdir}/scatterplot')

    def test_number_name(self):
        ScatterPlot(self.settings).main(
            sample_coordinate_df=read_tsv(f'{self.indir}/sample-coordinate-number-name.tsv'),
            x_column='PC1',
            y_column='PC2',
            hue_column='Group',
            colors=['#1f77b4', '#ff7f0e', '#2ca02c'],
            output_prefix=f'{self.outdir}/scatterplot')

    def test_many_groups(self):
        self.settings.for_publication = False
        ScatterPlot(self.settings).main(
            sample_coordinate_df=read_tsv(f'{self.indir}/sample-coordinate-many-groups.tsv'),
            x_column='PC1',
            y_column='PC2',
            hue_column='Group',
            colors=[
                'red', 'green', 'blue', 'salmon', 'cyan', 'purple', 'yellow', 'black', 'orange', 'pink',
                'brown', 'gray', 'magenta', 'olive', 'navy', 'teal', 'lime', 'aqua', 'maroon', 'fuchsia',
                'silver'
            ],
            output_prefix=f'{self.outdir}/scatterplot')


def read_tsv(tsv: str) -> pd.DataFrame:
    return pd.read_csv(tsv, sep='\t', index_col=0)
