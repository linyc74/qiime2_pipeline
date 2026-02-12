from os.path import exists
from qiime2_pipeline.picrust2 import PICRUSt2
from .setup import TestCase


class TestPICRUSt2(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = PICRUSt2(self.settings).main(
            labeled_feature_sequence_fa=f'{self.indir}/labeled-feature-sequence.fa',
            labeled_feature_table_tsv=f'{self.indir}/labeled-feature-table.tsv'
        )
        expected = {
            'picrust2-pathway': f'{self.outdir}/picrust2/picrust2-pathway-table.tsv',
            'picrust2-EC': f'{self.outdir}/picrust2/picrust2-EC-table.tsv',
            'picrust2-KEGG-ortholog': f'{self.outdir}/picrust2/picrust2-KEGG-ortholog-table.tsv',
        }
        self.assertDictEqual(expected, actual)
        for fpath in actual.values():
            with self.subTest(fpath=fpath):
                self.assertTrue(exists(fpath))

    def test_fungi_should_fail(self):
        actual = PICRUSt2(self.settings).main(
            labeled_feature_sequence_fa=f'{self.indir}/fungi-labeled-feature-sequence.fa',
            labeled_feature_table_tsv=f'{self.indir}/fungi-labeled-feature-table.tsv'
        )
        self.assertDictEqual({}, actual)

    def test_add_descriptions_pathway(self):
        PICRUSt2(self.settings).add_descriptions(
            in_tsv=f'{self.indir}/picrust2/pathways_out/path_abun_unstrat.tsv.gz',
            map_type='METACYC',
            out_tsv=f'{self.outdir}/picrust2-pathway-table.tsv'
        )

    def test_add_descriptions_ec(self):
        PICRUSt2(self.settings).add_descriptions(
            in_tsv=f'{self.indir}/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz',
            map_type='EC',
            out_tsv=f'{self.outdir}/picrust2-EC-table.tsv'
        )

    def test_add_descriptions_kegg_ortholog(self):
        PICRUSt2(self.settings).add_descriptions(
            in_tsv=f'{self.indir}/picrust2/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz',
            map_type='KO',
            out_tsv=f'{self.outdir}/picrust2-KEGG-ortholog-table.tsv'
        )


