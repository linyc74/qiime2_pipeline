from .setup import TestCase
from qiime2_pipeline.lefse_plot_cladogram import LefSePlotCladogram


class TestLefSePlotCladogram(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        LefSePlotCladogram().main(
            input_file=f'{self.indir}/lefse-genus-result.txt',
            output_file=f'{self.outdir}/lefse-genus-result.png',
            colormap='Set1',
        )
