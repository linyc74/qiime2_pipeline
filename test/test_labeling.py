from .setup import TestCase
from qiime2_pipeline.labeling import FeatureLabeling, add_genus_prefix_to_unidentified_species, \
    ensure_one_space_after_semicolon


class TestFeatureLabeling(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test_nb_bacteria(self):
        table_tsv, table_qza, sequence_fa, sequence_qza = FeatureLabeling(self.settings).main(
            taxonomy_qza=f'{self.indir}/taxonomy.qza',
            feature_table_qza=f'{self.indir}/dada2-feature-table.qza',
            feature_sequence_qza=f'{self.indir}/dada2-feature-sequence.qza',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            skip_otu=True)

        for expected, actual in [
            (f'{self.outdir}/labeled-feature-table.tsv', table_tsv),
            (f'{self.workdir}/labeled-feature-table.qza', table_qza),
            (f'{self.outdir}/labeled-feature-sequence.fa', sequence_fa),
            (f'{self.workdir}/labeled-feature-sequence.qza', sequence_qza),
        ]:
            self.assertFileExists(expected, actual)

    def test_nb_fungi(self):
        table_tsv, table_qza, sequence_fa, sequence_qza = FeatureLabeling(self.settings).main(
            taxonomy_qza=f'{self.indir}/fungi-taxonomy.qza',
            feature_table_qza=f'{self.indir}/fungi-dada2-feature-table.qza',
            feature_sequence_qza=f'{self.indir}/fungi-dada2-feature-sequence.qza',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            skip_otu=True)

        for expected, actual in [
            (f'{self.outdir}/labeled-feature-table.tsv', table_tsv),
            (f'{self.workdir}/labeled-feature-table.qza', table_qza),
            (f'{self.outdir}/labeled-feature-sequence.fa', sequence_fa),
            (f'{self.workdir}/labeled-feature-sequence.qza', sequence_qza),
        ]:
            self.assertFileExists(expected, actual)

    def test_vsearch_bacteria(self):
        table_tsv, table_qza, sequence_fa, sequence_qza = FeatureLabeling(self.settings).main(
            taxonomy_qza=f'{self.indir}/taxonomy-vsearch.qza',
            feature_table_qza=f'{self.indir}/dada2-feature-table.qza',
            feature_sequence_qza=f'{self.indir}/dada2-feature-sequence.qza',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            skip_otu=True)

        for expected, actual in [
            (f'{self.outdir}/labeled-feature-table.tsv', table_tsv),
            (f'{self.workdir}/labeled-feature-table.qza', table_qza),
            (f'{self.outdir}/labeled-feature-sequence.fa', sequence_fa),
            (f'{self.workdir}/labeled-feature-sequence.qza', sequence_qza),
        ]:
            self.assertFileExists(expected, actual)


class TestEnsureOneSpaceAfterSemicolon(TestCase):

    def test_already_with_space(self):
        actual = ensure_one_space_after_semicolon('k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella; s__intermedia')
        expected = 'k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella; s__intermedia'
        self.assertEqual(expected, actual)

    def test_no_space_after_semicolon(self):
        actual = ensure_one_space_after_semicolon('k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella;s__intermedia')
        expected = 'k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella; s__intermedia'
        self.assertEqual(expected, actual)


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
