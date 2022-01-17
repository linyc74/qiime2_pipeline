from .beta import BetaDiversity
from .alpha import AlphaDiversity
from .phylogeny import MafftFasttree
from .taxonomy import FeatureClassifier
from .template import Processor, Settings
from .generate_asv import FactoryGenerateASVCallable


class Qiime2Pipeline(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str
    nb_classifier_qza: str
    generate_asv_mode: str

    feature_sequence_qza: str
    feature_table_qza: str
    rooted_tree_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            nb_classifier_qza: str,
            generate_asv_mode: str):

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.nb_classifier_qza = nb_classifier_qza
        self.generate_asv_mode = generate_asv_mode

        self.generate_asv()
        self.taxonomic_classification()
        self.phylogenetic_tree()
        self.alpha_diversity()
        self.beta_diversity()

    def generate_asv(self):
        generate_asv = FactoryGenerateASVCallable(self.settings).main(
            mode=self.generate_asv_mode)

        self.feature_sequence_qza, self.feature_table_qza = generate_asv(
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def taxonomic_classification(self):
        FeatureClassifier(self.settings).main(
            representative_seq_qza=self.feature_sequence_qza,
            nb_classifier_qza=self.nb_classifier_qza)

    def phylogenetic_tree(self):
        self.rooted_tree_qza = MafftFasttree(self.settings).main(
            seq_qza=self.feature_sequence_qza)

    def alpha_diversity(self):
        AlphaDiversity(self.settings).main(
            feature_table_qza=self.feature_table_qza)

    def beta_diversity(self):
        BetaDiversity(self.settings).main(
            feature_table_qza=self.feature_table_qza,
            rooted_tree_qza=self.rooted_tree_qza)

