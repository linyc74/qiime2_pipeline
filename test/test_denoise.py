from .setup import TestCase
from qiime2_pipeline.denoise import Dada2SingleEnd, Dada2PairedEnd


class TestDada2SingleEnd(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        Dada2SingleEnd(self.settings).main(
            demultiplexed_seq_qza=f'{self.indir}/single-end-demultiplexed.qza',
            max_expected_error_bases=2.0,
        )


class TestDada2PairedEnd(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        Dada2PairedEnd(self.settings).main(
            demultiplexed_seq_qza=f'{self.indir}/paired-end-demultiplexed.qza',
            max_expected_error_bases=2.0,
        )
