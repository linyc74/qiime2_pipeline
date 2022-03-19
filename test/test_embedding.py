import pandas as pd
from os.path import exists
from qiime2_pipeline.embedding import PCA
from .setup import TestCase


class TestPCA(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        PCA(self.settings).main(
            feature_table_tsv=f'{self.indir}/feature-table.tsv'
        )
