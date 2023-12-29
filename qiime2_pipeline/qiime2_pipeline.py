from os import makedirs
from typing import List, Optional, Dict
from .lefse import LefSe
from .taxonomy import Taxonomy
from .template import Processor
from .phylogeny import Phylogeny
from .heatmap import PlotHeatmaps
from .alpha import AlphaDiversity
from .venn import PlotVennDiagrams
from .otu_clustering import Vsearch
from .taxon_table import TaxonTable
from .beta_my import MyBetaDiversity
from .labeling import FeatureLabeling
from .generate_asv import GenerateASV
from .beta_qiime import QiimeBetaDiversity
from .raw_read_counts import RawReadCounts
from .taxon_barplot import PlotTaxonBarplots
from .differential_abundance import DifferentialAbundance


class Qiime2Pipeline(Processor):

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: Optional[str]
    nb_classifier_qza: str
    paired_end_mode: str
    otu_identity: float
    skip_otu: bool
    classifier_reads_per_batch: int
    alpha_metrics: List[str]
    clip_r1_5_prime: int
    clip_r2_5_prime: int
    max_expected_error_bases: float
    heatmap_read_fraction: float
    n_taxa_barplot: int
    colormap: str

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
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: Optional[str],
            nb_classifier_qza: str,
            paired_end_mode: str,
            otu_identity: float,
            skip_otu: bool,
            classifier_reads_per_batch: int,
            alpha_metrics: List[str],
            clip_r1_5_prime: int,
            clip_r2_5_prime: int,
            max_expected_error_bases: float,
            heatmap_read_fraction: float,
            n_taxa_barplot: int,
            colormap: str):

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.nb_classifier_qza = nb_classifier_qza
        self.paired_end_mode = paired_end_mode
        self.otu_identity = otu_identity
        self.skip_otu = skip_otu
        self.classifier_reads_per_batch = classifier_reads_per_batch
        self.alpha_metrics = alpha_metrics
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime
        self.max_expected_error_bases = max_expected_error_bases
        self.heatmap_read_fraction = heatmap_read_fraction
        self.n_taxa_barplot = n_taxa_barplot
        self.colormap = colormap

        self.raw_read_counts()

        self.generate_asv()
        self.otu_clustering()

        self.taxonomic_classification()
        self.feature_labeling()

        self.phylogenetic_tree()

        self.alpha_diversity()

        self.qiime_beta_diversity()
        self.my_beta_diversity()

        self.taxon_table()
        self.plot_heatmaps()
        self.plot_venn_diagrams()
        self.taxon_barplot()
        self.lefse()
        self.differential_abundance()

        self.collect_log_files()

    def raw_read_counts(self):
        RawReadCounts(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def generate_asv(self):
        self.feature_table_qza, self.feature_sequence_qza = GenerateASV(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            paired_end_mode=self.paired_end_mode,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime,
            max_expected_error_bases=self.max_expected_error_bases)

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
            sample_sheet=self.sample_sheet,
            alpha_metrics=self.alpha_metrics,
            colormap=self.colormap)

    def qiime_beta_diversity(self):
        self.distance_matrix_tsvs = QiimeBetaDiversity(self.settings).main(
            feature_table_qza=self.labeled_feature_table_qza,
            rooted_tree_qza=self.rooted_tree_qza,
            sample_sheet=self.sample_sheet)

    def my_beta_diversity(self):
        MyBetaDiversity(self.settings).main(
            feature_table_tsv=self.labeled_feature_table_tsv,
            sample_sheet=self.sample_sheet)

    def taxon_table(self):
        self.taxon_table_tsv_dict = TaxonTable(self.settings).main(
            labeled_feature_table_tsv=self.labeled_feature_table_tsv)

    def plot_heatmaps(self):
        tsvs = [self.labeled_feature_table_tsv] + [v for v in self.taxon_table_tsv_dict.values()]
        PlotHeatmaps(self.settings).main(
            tsvs=tsvs,
            heatmap_read_fraction=self.heatmap_read_fraction)

    def plot_venn_diagrams(self):
        tsvs = [self.labeled_feature_table_tsv] + [v for v in self.taxon_table_tsv_dict.values()]
        PlotVennDiagrams(self.settings).main(
            tsvs=tsvs,
            sample_sheet=self.sample_sheet)

    def taxon_barplot(self):
        PlotTaxonBarplots(self.settings).main(
            taxon_table_tsv_dict=self.taxon_table_tsv_dict,
            n_taxa=self.n_taxa_barplot)

    def lefse(self):
        LefSe(self.settings).main(
            taxon_table_tsv_dict=self.taxon_table_tsv_dict,
            sample_sheet=self.sample_sheet,
            colormap=self.colormap)

    def differential_abundance(self):
        DifferentialAbundance(self.settings).main(
            taxon_table_tsv_dict=self.taxon_table_tsv_dict,
            sample_sheet=self.sample_sheet,
            colormap=self.colormap)

    def collect_log_files(self):
        makedirs(f'{self.outdir}/log', exist_ok=True)
        cmd = f'mv {self.outdir}/*.log {self.outdir}/log/'
        self.call(cmd)
