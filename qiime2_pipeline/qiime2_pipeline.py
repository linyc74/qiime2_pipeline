from typing import List, Optional, Dict
from .lefse import LefSe
from .taxonomy import Taxonomy
from .template import Processor
from .beta import BetaDiversity
from .phylogeny import Phylogeny
from .heatmap import PlotHeatmaps
from .alpha import AlphaDiversity
from .venn import PlotVennDiagrams
from .otu_clustering import Vsearch
from .taxon_table import TaxonTable
from .labeling import FeatureLabeling
from .taxon_barplot import PlotTaxonBarplots
from .feature_embedding import PCAProcess, NMDSProcess, TSNEProcess
from .generate_asv import FactoryGenerateASVPairedEnd, GenerateASVSingleEnd
from .beta_embedding import BatchPCoAProcess, BatchNMDSProcess, BatchTSNEProcess


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
    n_taxa_barplot: int

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
            n_taxa_barplot: int):

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
        self.n_taxa_barplot = n_taxa_barplot

        self.generate_asv()
        self.otu_clustering()
        self.taxonomic_classification()
        self.feature_labeling()
        self.phylogenetic_tree()
        self.alpha_diversity()
        self.beta_diversity()
        self.beta_embedding()
        self.taxon_table()
        self.lefse()
        self.plot_heatmaps()
        self.feature_embedding()
        self.taxon_barplot()
        self.plot_venn_diagrams()

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
        self.rooted_tree_qza = Phylogeny(self.settings).main(
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

    def beta_embedding(self):
        for Batch in [BatchPCoAProcess, BatchNMDSProcess, BatchTSNEProcess]:
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
        tsvs = [self.labeled_feature_table_tsv] + [v for v in self.taxon_table_tsv_dict.values()]
        PlotHeatmaps(self.settings).main(
            tsvs=tsvs,
            heatmap_read_fraction=self.heatmap_read_fraction)

    def feature_embedding(self):
        for Process in [PCAProcess, NMDSProcess, TSNEProcess]:
            Process(self.settings).main(
                tsv=self.labeled_feature_table_tsv,
                group_keywords=self.group_keywords)

    def taxon_barplot(self):
        PlotTaxonBarplots(self.settings).main(
            taxon_table_tsv_dict=self.taxon_table_tsv_dict,
            n_taxa=self.n_taxa_barplot)

    def plot_venn_diagrams(self):
        tsvs = [self.labeled_feature_table_tsv] + [v for v in self.taxon_table_tsv_dict.values()]
        PlotVennDiagrams(self.settings).main(
            tsvs=tsvs,
            group_keywords=self.group_keywords)
