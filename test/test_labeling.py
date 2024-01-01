from .setup import TestCase
from qiime2_pipeline.labeling import FeatureLabeling, add_genus_prefix_to_unidentified_species


class TestFeatureLabeling(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        table_tsv, table_qza, sequence_qza = FeatureLabeling(self.settings).main(
            taxonomy_qza=f'{self.indir}/taxonomy.qza',
            feature_table_qza=f'{self.indir}/dada2-feature-table.qza',
            feature_sequence_qza=f'{self.indir}/dada2-feature-sequence.qza',
            skip_otu=True)

        for expected, actual in [
            (f'{self.outdir}/labeled-feature-table.tsv', table_tsv),
            (f'{self.workdir}/labeled-feature-table.qza', table_qza),
            (f'{self.workdir}/labeled-feature-sequence.qza', sequence_qza),
        ]:
            self.assertFileExists(expected, actual)


class TestAddGenusPrefixToUnidentifiedSpecies(TestCase):

    def test_no_species_name(self):
        taxon = 'd__Bacteria'
        actual = add_genus_prefix_to_unidentified_species(taxon=taxon)
        self.assertEqual(taxon, actual)

    def test_proper_species_name(self):
        # species name contains genus name
        taxon = 'd__Bacteria; p__Fusobacteriota; c__Fusobacteriia; o__Fusobacteriales; f__Fusobacteriaceae; g__Fusobacterium; s__Fusobacterium_nucleatum'
        actual = add_genus_prefix_to_unidentified_species(taxon=taxon)
        self.assertEqual(taxon, actual)

    def test_species_name_without_genus(self):
        taxon = 'd__Bacteria; p__Fusobacteriota; c__Fusobacteriia; o__Fusobacteriales; f__Fusobacteriaceae; g__Fusobacterium; s__uncultured_bacterium'
        actual = add_genus_prefix_to_unidentified_species(taxon=taxon)
        expected = 'd__Bacteria; p__Fusobacteriota; c__Fusobacteriia; o__Fusobacteriales; f__Fusobacteriaceae; g__Fusobacterium; s__Fusobacterium_uncultured_bacterium'
        self.assertEqual(expected, actual)
