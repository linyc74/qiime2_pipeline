from os import makedirs
from typing import List, Optional, Dict
from .lefse import LefSe
from .tools import edit_fpath
from .taxonomy import Taxonomy
from .picrust2 import PICRUSt2
from .template import Processor
from .grouping import GetColors
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
    beta_diversity_feature_level: str
    colormap: str
    invert_colors: bool
    skip_differential_abundance: bool
    skip_picrust2: bool

    colors: list
    feature_table_qza: str
    feature_sequence_qza: str
    taxonomy_qza: str
    labeled_feature_table_tsv: str
    labeled_feature_table_qza: str
    labeled_feature_sequence_fa: str
    labeled_feature_sequence_qza: str
    rooted_tree_qza: str
    taxon_table_tsv_dict: Dict[str, str]
    picrust2_pathway_table_tsv: Optional[str]

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
            beta_diversity_feature_level: str,
            colormap: str,
            invert_colors: bool,
            skip_differential_abundance: bool,
            skip_picrust2: bool):

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
        self.beta_diversity_feature_level = beta_diversity_feature_level
        self.colormap = colormap
        self.invert_colors = invert_colors
        self.skip_differential_abundance = skip_differential_abundance
        self.skip_picrust2 = skip_picrust2

        self.copy_sample_sheet()

        self.raw_read_counts()
        self.set_colors()

        self.generate_asv()
        self.otu_clustering()

        self.taxonomic_classification()
        self.feature_labeling()
        self.taxon_table()
        self.picrust2()

        self.phylogenetic_tree()

        self.alpha_diversity()

        self.beta_diversity()

        self.plot_heatmaps()
        self.plot_venn_diagrams()
        self.taxon_barplot()
        self.lefse()
        self.differential_abundance()

        self.collect_log_files()

    def copy_sample_sheet(self):
        # to avoid the user from modifying or deleting the original sample sheet
        self.call(f'cp "{self.sample_sheet}" "{self.workdir}/"')
        self.sample_sheet = edit_fpath(
            fpath=self.sample_sheet,
            dstdir=self.workdir)

    def raw_read_counts(self):
        RawReadCounts(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def set_colors(self):
        self.colors = GetColors(self.settings).main(
            sample_sheet=self.sample_sheet,
            colormap=self.colormap,
            invert_colors=self.invert_colors)

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
        self.labeled_feature_table_tsv, self.labeled_feature_table_qza, \
            self.labeled_feature_sequence_fa, self.labeled_feature_sequence_qza = FeatureLabeling(self.settings).main(
                taxonomy_qza=self.taxonomy_qza,
                feature_table_qza=self.feature_table_qza,
                feature_sequence_qza=self.feature_sequence_qza,
                sample_sheet=self.sample_sheet,
                skip_otu=self.skip_otu)

    def taxon_table(self):
        self.taxon_table_tsv_dict = TaxonTable(self.settings).main(
            labeled_feature_table_tsv=self.labeled_feature_table_tsv)

    def picrust2(self):
        if self.skip_picrust2:
            self.picrust2_pathway_table_tsv = None
            return
        self.picrust2_pathway_table_tsv = PICRUSt2(self.settings).main(
            labeled_feature_sequence_fa=self.labeled_feature_sequence_fa,
            labeled_feature_table_tsv=self.labeled_feature_table_tsv)

    def phylogenetic_tree(self):
        self.rooted_tree_qza = Phylogeny(self.settings).main(
            seq_qza=self.labeled_feature_sequence_qza)

    def alpha_diversity(self):
        AlphaDiversity(self.settings).main(
            feature_table_qza=self.labeled_feature_table_qza,
            sample_sheet=self.sample_sheet,
            alpha_metrics=self.alpha_metrics,
            colors=self.colors)

    def beta_diversity(self):
        feature_table_tsv = self.labeled_feature_table_tsv if self.beta_diversity_feature_level == 'feature' \
            else self.taxon_table_tsv_dict[self.beta_diversity_feature_level]

        QiimeBetaDiversity(self.settings).main(
            feature_table_tsv=feature_table_tsv,
            rooted_tree_qza=self.rooted_tree_qza,
            sample_sheet=self.sample_sheet,
            colors=self.colors)

        MyBetaDiversity(self.settings).main(
            feature_table_tsv=feature_table_tsv,
            sample_sheet=self.sample_sheet,
            colors=self.colors)

    def plot_heatmaps(self):
        tsvs = [self.labeled_feature_table_tsv] + [v for v in self.taxon_table_tsv_dict.values()]
        PlotHeatmaps(self.settings).main(
            tsvs=tsvs,
            heatmap_read_fraction=self.heatmap_read_fraction,
            sample_sheet=self.sample_sheet)

    def plot_venn_diagrams(self):
        tsvs = [self.labeled_feature_table_tsv] + [v for v in self.taxon_table_tsv_dict.values()]
        PlotVennDiagrams(self.settings).main(
            tsvs=tsvs,
            sample_sheet=self.sample_sheet,
            colors=self.colors)

    def taxon_barplot(self):
        PlotTaxonBarplots(self.settings).main(
            taxon_table_tsv_dict=self.taxon_table_tsv_dict,
            n_taxa=self.n_taxa_barplot,
            sample_sheet=self.sample_sheet)

    def lefse(self):
        table_tsv_dict = self.taxon_table_tsv_dict.copy()
        if self.picrust2_pathway_table_tsv is not None:
            table_tsv_dict['picrust2-pathway'] = self.picrust2_pathway_table_tsv
        LefSe(self.settings).main(
            table_tsv_dict=table_tsv_dict,
            sample_sheet=self.sample_sheet,
            colors=self.colors)

    def differential_abundance(self):
        if self.skip_differential_abundance:
            return
        DifferentialAbundance(self.settings).main(
            taxon_table_tsv_dict=self.taxon_table_tsv_dict,
            sample_sheet=self.sample_sheet,
            colors=self.colors)

    def collect_log_files(self):
        makedirs(f'{self.outdir}/log', exist_ok=True)
        cmd = f'mv "{self.outdir}"/*.log "{self.outdir}"/log/'
        self.call(cmd)
