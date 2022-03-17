from typing import List, Optional, Dict
from .lefse import LefSe
from .taxonomy import Taxonomy
from .template import Processor
from .beta import BetaDiversity
from .heatmap import PlotHeatmaps
from .alpha import AlphaDiversity
from .otu_clustering import Vsearch
from .taxon_table import TaxonTable
from .phylogeny import MafftFasttree
from .labeling import FeatureLabeling
from .normalization import CountNormalization
from .beta_embedding import BatchPCoA, BatchNMDS, BatchTSNE
from .generate_asv import FactoryGenerateASVPairedEnd, GenerateASVSingleEnd


class Qiime2Pipeline(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: Optional[str]
    nb_classifier_qza: str
    paired_end_mode: str
    group_keywords: List[str]
    otu_identity: float
    skip_otu: bool
    classifier_reads_per_batch: int
    alpha_metrics: List[str]
    clip_r1_5_prime: int
    clip_r2_5_prime: int
    heatmap_read_fraction: float
    log_pseudocount: bool

    feature_table_qza: str
    feature_sequence_qza: str
    taxonomy_qza: str
    labeled_feature_table_tsv: str
    labeled_feature_table_qza: str
    labeled_feature_sequence_qza: str
    rooted_tree_qza: str
    distance_matrix_tsvs: List[str]
    taxon_table_tsv_dict: Dict[str, str]

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: Optional[str],
            nb_classifier_qza: str,
            paired_end_mode: str,
            group_keywords: List[str],
            otu_identity: float,
            skip_otu: bool,
            classifier_reads_per_batch: int,
            alpha_metrics: List[str],
            clip_r1_5_prime: int,
            clip_r2_5_prime: int,
            heatmap_read_fraction: float,
            log_pseudocount: bool):

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.nb_classifier_qza = nb_classifier_qza
        self.paired_end_mode = paired_end_mode
        self.group_keywords = group_keywords
        self.otu_identity = otu_identity
        self.skip_otu = skip_otu
        self.classifier_reads_per_batch = classifier_reads_per_batch
        self.alpha_metrics = alpha_metrics
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime
        self.heatmap_read_fraction = heatmap_read_fraction
        self.log_pseudocount = log_pseudocount

        self.generate_asv()
        self.otu_clustering()
        self.taxonomic_classification()
        self.feature_labeling()
        self.phylogenetic_tree()
        self.alpha_diversity()
        self.beta_diversity()
        self.dimensionality_reduction()
        self.taxon_table()
        self.lefse()
        self.plot_heatmaps()

    def generate_asv(self):
        if self.fq2_suffix is None:
            self.feature_table_qza, self.feature_sequence_qza = GenerateASVSingleEnd(self.settings).main(
                fq_dir=self.fq_dir,
                fq_suffix=self.fq1_suffix,
                clip_5_prime=self.clip_r1_5_prime)
        else:
            generate_asv = FactoryGenerateASVPairedEnd(self.settings).main(
                paired_end_mode=self.paired_end_mode)
            self.feature_table_qza, self.feature_sequence_qza = generate_asv(
                fq_dir=self.fq_dir,
                fq1_suffix=self.fq1_suffix,
                fq2_suffix=self.fq2_suffix,
                clip_r1_5_prime=self.clip_r1_5_prime,
                clip_r2_5_prime=self.clip_r2_5_prime)

    def otu_clustering(self):
        if self.skip_otu:
            return
        self.feature_table_qza, self.feature_sequence_qza = Vsearch(self.settings).main(
            feature_table_qza=self.feature_table_qza,
            feature_sequence_qza=self.feature_sequence_qza,
            identity=self.otu_identity)

    def taxonomic_classification(self):
        self.taxonomy_qza = Taxonomy(self.settings).main(
            representative_seq_qza=self.feature_sequence_qza,
            nb_classifier_qza=self.nb_classifier_qza,
            classifier_reads_per_batch=self.classifier_reads_per_batch)

    def feature_labeling(self):
        self.labeled_feature_table_tsv, \
            self.labeled_feature_table_qza, \
            self.labeled_feature_sequence_qza = FeatureLabeling(self.settings).main(
                taxonomy_qza=self.taxonomy_qza,
                feature_table_qza=self.feature_table_qza,
                feature_sequence_qza=self.feature_sequence_qza,
                skip_otu=self.skip_otu)

    def phylogenetic_tree(self):
        self.rooted_tree_qza = MafftFasttree(self.settings).main(
            seq_qza=self.labeled_feature_sequence_qza)

    def alpha_diversity(self):
        AlphaDiversity(self.settings).main(
            feature_table_qza=self.labeled_feature_table_qza,
            group_keywords=self.group_keywords,
            alpha_metrics=self.alpha_metrics)

    def beta_diversity(self):
        self.distance_matrix_tsvs = BetaDiversity(self.settings).main(
            feature_table_qza=self.labeled_feature_table_qza,
            rooted_tree_qza=self.rooted_tree_qza)

    def dimensionality_reduction(self):
        for Batch in [BatchPCoA, BatchNMDS, BatchTSNE]:
            Batch(self.settings).main(
                distance_matrix_tsvs=self.distance_matrix_tsvs,
                group_keywords=self.group_keywords)

    def taxon_table(self):
        self.taxon_table_tsv_dict = TaxonTable(self.settings).main(
            labeled_feature_table_tsv=self.labeled_feature_table_tsv)

    def lefse(self):
        LefSe(self.settings).main(
            taxon_table_tsv_dict=self.taxon_table_tsv_dict,
            group_keywords=self.group_keywords)

    def plot_heatmaps(self):
        PlotHeatmaps(self.settings).main(
            tsvs=list(self.taxon_table_tsv_dict.values()),
            heatmap_read_fraction=self.heatmap_read_fraction)
