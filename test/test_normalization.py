from .setup import TestCase
from qiime2_pipeline.normalization import FeatureNormalization


class TestFeatureNormalization(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        table_tsv, table_qza = FeatureNormalization(self.settings).main(
            feature_table_qza=f'{self.indir}/table.qza',
            log_pseudocount=True
        )

        for expected, actual in [
            (f'{self.outdir}/table-normalized.tsv', table_tsv),
            (f'{self.workdir}/table-normalized.qza', table_qza),
        ]:
            self.assertFileExists(expected, actual)
