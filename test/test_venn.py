from qiime2_pipeline.venn import PlotVennDiagrams, ProcessTsvPlotVenn
from .setup import TestCase


class TestPlotVennDiagrams(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        tsvs = [
            f'{self.indir}/tsvs/phylum-table.tsv',
            f'{self.indir}/tsvs/class-table.tsv',
            f'{self.indir}/tsvs/order-table.tsv',
            f'{self.indir}/tsvs/family-table.tsv',
            f'{self.indir}/tsvs/genus-table.tsv',
            f'{self.indir}/tsvs/species-table.tsv',
        ]
        PlotVennDiagrams(self.settings).main(
            tsvs=tsvs,
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)]
        )


class TestProcessTsvPlotVenn(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_1_group(self):
        with self.assertRaises(AssertionError):
            ProcessTsvPlotVenn(self.settings).main(
                tsv=f'{self.indir}/mock-feature-table.tsv',
                sample_sheet=f'{self.indir}/mock-sample-sheet-1-group.csv',
                colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)],
                dstdir=self.outdir)

    def test_2_groups(self):
        ProcessTsvPlotVenn(self.settings).main(
            tsv=f'{self.indir}/mock-feature-table.tsv',
            sample_sheet=f'{self.indir}/mock-sample-sheet-2-groups.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)],
            dstdir=self.outdir)

    def test_3_groups(self):
        ProcessTsvPlotVenn(self.settings).main(
            tsv=f'{self.indir}/mock-feature-table.tsv',
            sample_sheet=f'{self.indir}/mock-sample-sheet-3-groups.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0), (0.3, 0.6, 0.0, 1.0)],
            dstdir=self.outdir)

    def test_4_groups(self):
        with self.assertRaises(AssertionError):
            ProcessTsvPlotVenn(self.settings).main(
                tsv=f'{self.indir}/mock-feature-table.tsv',
                sample_sheet=f'{self.indir}/mock-sample-sheet-4-groups.csv',
                colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)],
                dstdir=self.outdir)
