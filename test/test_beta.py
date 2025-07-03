import pandas as pd
from os.path import exists
from .setup import TestCase
from qiime2_pipeline.beta import BetaDiversity, RunAllBetaMetricsToDistanceMatrixTsvs, PCoAProcess, PCAProcess, PCACore, ScatterPlot


class TestQiimeBetaDiversity(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        BetaDiversity(self.settings).main(
            feature_table_tsv=f'{self.indir}/feature-table.tsv',
            rooted_tree_qza=f'{self.indir}/fasttree-rooted.qza',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0),  'green'],
        )


class TestRunAllBetaMetricsToTsvs(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = RunAllBetaMetricsToDistanceMatrixTsvs(self.settings).main(
            feature_table_qza=f'{self.indir}/feature-table.qza',
            rooted_tree_qza=f'{self.indir}/fasttree-rooted.qza'
        )
        expected = [
            f'{self.outdir}/beta-diversity/jaccard.tsv',
            f'{self.outdir}/beta-diversity/euclidean.tsv',
            f'{self.outdir}/beta-diversity/braycurtis.tsv',
            f'{self.outdir}/beta-diversity/weighted_unifrac.tsv',
            f'{self.outdir}/beta-diversity/weighted_normalized_unifrac.tsv',
            f'{self.outdir}/beta-diversity/generalized_unifrac.tsv',
            f'{self.outdir}/beta-diversity/unweighted_unifrac.tsv',
        ]
        for e, a in zip(expected, actual):
            self.assertFileExists(expected=e, actual=a)


class TestPCoAProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        name = 'jaccard'
        PCoAProcess(self.settings).main(
            tsv=f'{self.indir}/distance-matrix-tsvs/{name}.tsv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0),  'green'],
        )
        for f in [
            f'{self.outdir}/beta-embedding/{name}-pcoa-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/{name}-pcoa-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/{name}-pcoa-sample-coordinate.pdf',
            f'{self.outdir}/beta-embedding/{name}-pcoa-proportion-explained.tsv',
        ]:
            self.assertTrue(exists(f))


class TestPCAProcess(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        PCAProcess(self.settings).main(
            tsv=f'{self.indir}/feature-table-large.tsv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0), (0.1, 0.9, 0.3, 1.0)],
        )
        for f in [
            f'{self.outdir}/beta-embedding/feature-table-large-pca-sample-coordinate.tsv',
            f'{self.outdir}/beta-embedding/feature-table-large-pca-sample-coordinate.png',
            f'{self.outdir}/beta-embedding/feature-table-large-pca-sample-coordinate.pdf',
            f'{self.outdir}/beta-embedding/feature-table-large-pca-proportion-explained.tsv',
        ]:
            self.assertTrue(exists(f))


class TestPCACore(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        sample_coordinate_df, proportion_explained_series = PCACore(self.settings).main(
            df=read_tsv(f'{self.indir}/feature-table-large.tsv'),
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
            x_label_suffix=' (0.0%)',
            y_label_suffix=' (0.0%)',
            output_prefix=f'{self.outdir}/scatterplot')

    def test_long_name(self):
        ScatterPlot(self.settings).main(
            sample_coordinate_df=read_tsv(f'{self.indir}/sample-coordinate-long-name.tsv'),
            x_column='PC1',
            y_column='PC2',
            hue_column='Group',
            colors=['#1f77b4', '#ff7f0e', '#2ca02c'],
            x_label_suffix=' (0.0%)',
            y_label_suffix=' (0.0%)',
            output_prefix=f'{self.outdir}/scatterplot')

    def test_number_name(self):
        ScatterPlot(self.settings).main(
            sample_coordinate_df=read_tsv(f'{self.indir}/sample-coordinate-number-name.tsv'),
            x_column='PC1',
            y_column='PC2',
            hue_column='Group',
            colors=['#1f77b4', '#ff7f0e', '#2ca02c'],
            x_label_suffix=' (0.0%)',
            y_label_suffix=' (0.0%)',
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
            x_label_suffix=' (0.0%)',
            y_label_suffix=' (0.0%)',
            output_prefix=f'{self.outdir}/scatterplot')


def read_tsv(tsv: str) -> pd.DataFrame:
    return pd.read_csv(tsv, sep='\t', index_col=0)
