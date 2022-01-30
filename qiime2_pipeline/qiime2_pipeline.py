from typing import List
from .taxonomy import Taxonomy
from .beta import BetaDiversity
from .alpha import AlphaDiversity
from .phylogeny import MafftFasttree
from .labeling import FeatureLabeling
from .template import Processor, Settings
from .generate_asv import FactoryGenerateASVCallable
from .dim_reduction import BatchPCoA, BatchNMDS, BatchTSNE


class Qiime2Pipeline(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str
    nb_classifier_qza: str
    paired_end_mode: str
    group_keywords: List[str]

    feature_sequence_qza: str
    feature_table_qza: str
    taxonomy_qza: str
    labeled_feature_sequence_qza: str
    labeled_feature_table_qza: str
    rooted_tree_qza: str
    distance_matrix_tsvs: List[str]

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            nb_classifier_qza: str,
            paired_end_mode: str,
            group_keywords: List[str]):

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.nb_classifier_qza = nb_classifier_qza
        self.paired_end_mode = paired_end_mode
        self.group_keywords = group_keywords

        self.generate_asv()
        self.taxonomic_classification()
        self.feature_labeling()
        self.phylogenetic_tree()
        self.alpha_diversity()
        self.beta_diversity()
        self.dimensionality_reduction()

    def generate_asv(self):
        generate_asv = FactoryGenerateASVCallable(self.settings).main(
            paired_end_mode=self.paired_end_mode)

        self.feature_table_qza, self.feature_sequence_qza = generate_asv(
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def taxonomic_classification(self):
        self.taxonomy_qza = Taxonomy(self.settings).main(
            representative_seq_qza=self.feature_sequence_qza,
            nb_classifier_qza=self.nb_classifier_qza)

    def feature_labeling(self):
        self.labeled_feature_table_qza, self.labeled_feature_sequence_qza \
            = FeatureLabeling(self.settings).main(
                taxonomy_qza=self.taxonomy_qza,
                feature_table_qza=self.feature_table_qza,
                feature_sequence_qza=self.feature_sequence_qza)

    def phylogenetic_tree(self):
        self.rooted_tree_qza = MafftFasttree(self.settings).main(
            seq_qza=self.labeled_feature_sequence_qza)

    def alpha_diversity(self):
        AlphaDiversity(self.settings).main(
            feature_table_qza=self.labeled_feature_table_qza)

    def beta_diversity(self):
        self.distance_matrix_tsvs = BetaDiversity(self.settings).main(
            feature_table_qza=self.labeled_feature_table_qza,
            rooted_tree_qza=self.rooted_tree_qza)

    def dimensionality_reduction(self):
        for Batch in [BatchPCoA, BatchNMDS, BatchTSNE]:
            Batch(self.settings).main(self.distance_matrix_tsvs)
