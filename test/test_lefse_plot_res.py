from .setup import TestCase
from qiime2_pipeline.lefse_plot_res import LefSePlotRes


class TestLefSePlotRes(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        LefSePlotRes().main(
            input_file=f'{self.indir}/lefse-genus-result.txt',
            output_file=f'{self.outdir}/lefse-genus-result.png',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )
