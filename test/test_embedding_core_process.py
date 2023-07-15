import pandas as pd
from qiime2_pipeline.embedding_core_process import PCACore, NMDSCore, TSNECore
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

    def test_nmds_core(self):
        sample_coordinate_df, stress = NMDSCore(self.settings).main(
            df=self.feature_table_df,
            data_structure='row_features'
        )
        expected = read_tsv(f'{self.indir}/nmds_sample_coordinate_df.tsv')
        self.assertDataFrameEqual(expected, sample_coordinate_df)

    def test_tsne_core(self):
        sample_coordinate_df = TSNECore(self.settings).main(
            df=self.feature_table_df,
            data_structure='row_features'
        )

        # save tsv to avoid floating point imprecision error
        actual = f'{self.outdir}/tsne_sample_coordinate_df.tsv'
        expected = f'{self.indir}/tsne_sample_coordinate_df.tsv'
        sample_coordinate_df.to_csv(actual, sep='\t')
        self.assertFileEqual(expected, actual)


def read_tsv(tsv: str) -> pd.DataFrame:
    return pd.read_csv(tsv, sep='\t', index_col=0)
