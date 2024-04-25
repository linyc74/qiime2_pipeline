import shutil
from qiime2_pipeline.picrust2 import PICRUSt2
from .setup import TestCase


class TestPICRUSt2(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

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
        for e, a in zip(expected, actual):
            self.assertFileExists(e, a)

    def test_add_description_pathway(self):
        PICRUSt2(self.settings).add_description(
            in_tsv=f'{self.indir}/picrust2/pathways_out/path_abun_unstrat.tsv.gz',
            map_type='METACYC',
            out_tsv=f'{self.outdir}/picrust2-pathway-table.tsv'
        )

    def test_add_description_ec(self):
        PICRUSt2(self.settings).add_description(
            in_tsv=f'{self.indir}/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz',
            map_type='EC',
            out_tsv=f'{self.outdir}/picrust2-EC-table.tsv'
        )

    def test_add_description_kegg_ortholog(self):
        PICRUSt2(self.settings).add_description(
            in_tsv=f'{self.indir}/picrust2/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz',
            map_type='KO',
            out_tsv=f'{self.outdir}/picrust2-KEGG-ortholog-table.tsv'
        )


