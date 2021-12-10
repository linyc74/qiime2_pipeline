from os.path import basename
from .beta import BetaDiversity
from .denoise import Dada2Paired
from .alpha import AlphaDiversity
from .phylogeny import MafftFasttree
from .taxonomy import FeatureClassifier
from .trimming import CutadaptTrimPaired
from .template import Processor, Settings
from .importing import ImportPairedEndFastq
from .exporting import ExportTaxonomy, ExportTree, ExportFeatureTable, ExportFeatureSequence, ExportAlignedSequence


class Qiime2Pipeline(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str
    nb_classifier_qza: str
    read_length: int

    untrimmed_reads_qza: str
    trimmed_reads_qza: str
    representative_seq_qza: str
    feature_table_qza: str
    taxonomy_qza: str
    aligned_seq_qza: str
    masked_aligned_seq_qza: str
    unrooted_tree_qza: str
    rooted_tree_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            nb_classifier_qza: str,
            read_length: int):

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.nb_classifier_qza = nb_classifier_qza
        self.read_length = read_length

        self.import_paired_end_fastq()
        self.trimming()
        self.denoise()
        self.taxonomic_classification()
        self.phylogenetic_tree()
        self.alpha_diversity()
        self.beta_diversity()
        self.export_data()

    def import_paired_end_fastq(self):
        self.untrimmed_reads_qza = ImportPairedEndFastq(self.settings).main(
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def trimming(self):
        self.trimmed_reads_qza = CutadaptTrimPaired(self.settings).main(
            untrimmed_reads_qza=self.untrimmed_reads_qza)

    def denoise(self):
        self.representative_seq_qza, self.feature_table_qza = Dada2Paired(self.settings).main(
            demultiplexed_seq_qza=self.trimmed_reads_qza,
            truncate_length=self.read_length)

    def taxonomic_classification(self):
        self.taxonomy_qza = FeatureClassifier(self.settings).main(
            representative_seq_qza=self.representative_seq_qza,
            nb_classifier_qza=self.nb_classifier_qza)

    def phylogenetic_tree(self):
        self.aligned_seq_qza, \
            self.masked_aligned_seq_qza, \
            self.unrooted_tree_qza, \
            self.rooted_tree_qza = MafftFasttree(self.settings).main(
                seq_qza=self.representative_seq_qza)

    def alpha_diversity(self):
        AlphaDiversity(self.settings).main(
            feature_table_qza=self.feature_table_qza)

    def beta_diversity(self):
        BetaDiversity(self.settings).main(
            feature_table_qza=self.feature_table_qza,
            rooted_tree_qza=self.rooted_tree_qza)

    def export_data(self):
        ExportData(self.settings).main(
            representative_seq_qza=self.representative_seq_qza,
            feature_table_qza=self.feature_table_qza,
            taxonomy_qza=self.taxonomy_qza,
            aligned_seq_qza=self.aligned_seq_qza,
            masked_aligned_seq_qza=self.masked_aligned_seq_qza,
            unrooted_tree_qza=self.unrooted_tree_qza,
            rooted_tree_qza=self.rooted_tree_qza)


class ExportData(Processor):

    representative_seq_qza: str
    feature_table_qza: str
    taxonomy_qza: str
    aligned_seq_qza: str
    masked_aligned_seq_qza: str
    unrooted_tree_qza: str
    rooted_tree_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            representative_seq_qza: str,
            feature_table_qza: str,
            taxonomy_qza: str,
            aligned_seq_qza: str,
            masked_aligned_seq_qza: str,
            unrooted_tree_qza: str,
            rooted_tree_qza: str):

        self.representative_seq_qza = representative_seq_qza
        self.feature_table_qza = feature_table_qza
        self.taxonomy_qza = taxonomy_qza
        self.aligned_seq_qza = aligned_seq_qza
        self.masked_aligned_seq_qza = masked_aligned_seq_qza
        self.unrooted_tree_qza = unrooted_tree_qza
        self.rooted_tree_qza = rooted_tree_qza

        self.export_feature_sequence()
        self.export_feature_table()
        self.export_taxonomy()
        self.export_aligned_sequences()
        self.export_trees()

    def export_feature_sequence(self):
        fname = basename(self.representative_seq_qza)[:-len('.qza')] + '.fa'
        ExportFeatureSequence(self.settings).main(
            feature_sequence_qza=self.representative_seq_qza,
            output_fa=f'{self.outdir}/{fname}'
        )

    def export_feature_table(self):
        fname = basename(self.feature_table_qza)[:-len('.qza')] + '.tsv'
        ExportFeatureTable(self.settings).main(
            feature_table_qza=self.feature_table_qza,
            output_tsv=f'{self.outdir}/{fname}'
        )

    def export_taxonomy(self):
        fname = basename(self.taxonomy_qza)[:-len('.qza')] + '.tsv'
        ExportTaxonomy(self.settings).main(
            taxonomy_qza=self.taxonomy_qza,
            output_tsv=f'{self.outdir}/{fname}'
        )

    def export_aligned_sequences(self):
        for qza in [self.aligned_seq_qza, self.masked_aligned_seq_qza]:
            fname = basename(qza)[:-len('.qza')] + '.fa'
            ExportAlignedSequence(self.settings).main(
                aligned_sequence_qza=qza,
                output_fa=f'{self.outdir}/{fname}'
            )

    def export_trees(self):
        for tree_qza in [self.rooted_tree_qza, self.unrooted_tree_qza]:
            fname = basename(tree_qza)[:-len('.qza')] + '.nwk'
            ExportTree(self.settings).main(
                tree_qza=tree_qza,
                output_nwk=f'{self.outdir}/{fname}'
            )
