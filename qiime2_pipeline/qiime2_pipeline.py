from os import makedirs
from typing import List, Optional, Dict
from .lefse import LefSe
from .taxonomy import Taxonomy
from .picrust2 import PICRUSt2
from .template import Processor
from .grouping import GetColors
from .phylogeny import Phylogeny
from .heatmap import PlotHeatmaps
from .alpha import AlphaDiversity
from .venn import PlotVennDiagrams
from .taxon_table import TaxonTable
from .beta_my import MyBetaDiversity
from .labeling import FeatureLabeling
from .generate_asv import GenerateASV
from .exporting import ExportFeatureTable
from .beta_qiime import QiimeBetaDiversity
from .raw_read_counts import RawReadCounts
from .taxon_barplot import PlotTaxonBarplots
from .sample_sheet import TranscribeSampleSheet
from .alpha_rarefaction import AlphaRarefaction
from .differential_abundance import DifferentialAbundance
from .generate_otu import GenerateOTU, GenerateNanoporeOTU


class Qiime2Pipeline(Processor):

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: Optional[str]

    sequencing_platform: str

    clip_r1_5_prime: int
    clip_r2_5_prime: int

    paired_end_mode: str
    max_expected_error_bases: float

    otu_identity: float
    skip_otu: bool

    feature_classifier: str
    nb_classifier_qza: Optional[str]
    classifier_reads_per_batch: int
    reference_sequence_qza: Optional[str]
    reference_taxonomy_qza: Optional[str]
    vsearch_classifier_max_hits: int

    alpha_metrics: List[str]
    beta_diversity_feature_level: str
    heatmap_read_fraction: float
    n_taxa_barplot: int
    colormap: str
    invert_colors: bool
    skip_differential_abundance: bool
    differential_abundance_p_value: float
    run_picrust2: bool

    colors: list
    feature_table_qza: str
    feature_sequence_qza: str
    taxonomy_qza: str
    labeled_feature_table_tsv: str
    labeled_feature_table_qza: str
    labeled_feature_sequence_fa: str
    labeled_feature_sequence_qza: str
    taxon_table_tsv_dict: Dict[str, str]
    picrust2_table_tsv_dict: Dict[str, str]

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: Optional[str],

            sequencing_platform: str,

            clip_r1_5_prime: int,
            clip_r2_5_prime: int,

            paired_end_mode: str,
            max_expected_error_bases: float,

            otu_identity: float,
            skip_otu: bool,

            feature_classifier: str,
            nb_classifier_qza: Optional[str],
            classifier_reads_per_batch: int,
            reference_sequence_qza: Optional[str],
            reference_taxonomy_qza: Optional[str],
            vsearch_classifier_max_hits: int,

            alpha_metrics: List[str],
            beta_diversity_feature_level: str,
            heatmap_read_fraction: float,
            n_taxa_barplot: int,
            colormap: str,
            invert_colors: bool,
            skip_differential_abundance: bool,
            differential_abundance_p_value: float,
            run_picrust2: bool):

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix

        self.sequencing_platform = sequencing_platform

        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        self.paired_end_mode = paired_end_mode
        self.max_expected_error_bases = max_expected_error_bases

        self.otu_identity = otu_identity
        self.skip_otu = skip_otu

        self.feature_classifier = feature_classifier
        self.nb_classifier_qza = nb_classifier_qza
        self.classifier_reads_per_batch = classifier_reads_per_batch
        self.reference_sequence_qza = reference_sequence_qza
        self.reference_taxonomy_qza = reference_taxonomy_qza
        self.vsearch_classifier_max_hits = vsearch_classifier_max_hits

        self.alpha_metrics = alpha_metrics
        self.beta_diversity_feature_level = beta_diversity_feature_level
        self.heatmap_read_fraction = heatmap_read_fraction
        self.n_taxa_barplot = n_taxa_barplot
        self.colormap = colormap
        self.invert_colors = invert_colors
        self.skip_differential_abundance = skip_differential_abundance
        self.differential_abundance_p_value = differential_abundance_p_value
        self.run_picrust2 = run_picrust2

        self.transcribe_sample_sheet()

        self.raw_read_counts()
        self.set_colors()

        self.generate_asv_otu()

        self.taxonomic_classification()
        self.feature_labeling()
        self.taxon_table()
        self.picrust2()

        self.alpha_diversity()
        self.alpha_rarefaction()

        self.phylogeny_and_beta_diversity()

        self.plot_heatmaps()
        self.plot_venn_diagrams()
        self.taxon_barplot()
        self.lefse()
        self.differential_abundance()

        self.collect_log_files()

    def transcribe_sample_sheet(self):
        self.sample_sheet = TranscribeSampleSheet(self.settings).main(
            sample_sheet=self.sample_sheet)

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

    def generate_asv_otu(self):
        if self.sequencing_platform == 'nanopore':
            self.feature_table_qza, self.feature_sequence_qza = GenerateNanoporeOTU(self.settings).main(
                sample_sheet=self.sample_sheet,
                fq_dir=self.fq_dir,
                fq_suffix=self.fq1_suffix,
                identity=self.otu_identity)

        elif self.sequencing_platform in ['illumina', 'pacbio']:
            self.feature_table_qza, self.feature_sequence_qza = GenerateASV(self.settings).main(
                sample_sheet=self.sample_sheet,
                fq_dir=self.fq_dir,
                fq1_suffix=self.fq1_suffix,
                fq2_suffix=self.fq2_suffix,
                pacbio=self.sequencing_platform == 'pacbio',
                paired_end_mode=self.paired_end_mode,
                clip_r1_5_prime=self.clip_r1_5_prime,
                clip_r2_5_prime=self.clip_r2_5_prime,
                max_expected_error_bases=self.max_expected_error_bases)
            if not self.skip_otu:
                self.feature_table_qza, self.feature_sequence_qza = GenerateOTU(self.settings).main(
                    feature_table_qza=self.feature_table_qza,
                    feature_sequence_qza=self.feature_sequence_qza,
                    identity=self.otu_identity)

        else:
            raise ValueError(f'Invalid sequencing platform: {self.sequencing_platform}')

    def taxonomic_classification(self):
        self.taxonomy_qza = Taxonomy(self.settings).main(
            representative_seq_qza=self.feature_sequence_qza,
            feature_classifier=self.feature_classifier,
            nb_classifier_qza=self.nb_classifier_qza,
            classifier_reads_per_batch=self.classifier_reads_per_batch,
            reference_sequence_qza=self.reference_sequence_qza,
            reference_taxonomy_qza=self.reference_taxonomy_qza,
            vsearch_classifier_max_hits=self.vsearch_classifier_max_hits)

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
        if self.run_picrust2:
            self.picrust2_table_tsv_dict = PICRUSt2(self.settings).main(
                labeled_feature_sequence_fa=self.labeled_feature_sequence_fa,
                labeled_feature_table_tsv=self.labeled_feature_table_tsv)
        else:
            self.picrust2_table_tsv_dict = {}

    def alpha_diversity(self):
        AlphaDiversity(self.settings).main(
            feature_table_qza=self.feature_table_qza,  # no need to use taxonomy-labeled feature table
            sample_sheet=self.sample_sheet,
            alpha_metrics=self.alpha_metrics,
            colors=self.colors)

    def alpha_rarefaction(self):
        AlphaRarefaction(self.settings).main(feature_table_qza=self.feature_table_qza)

    def phylogeny_and_beta_diversity(self):
        if self.beta_diversity_feature_level == 'feature':
            feature_table_tsv = ExportFeatureTable(self.settings).main(
                feature_table_qza=self.feature_table_qza)  # no need to use taxonomy-labeled feature table
        else:
            feature_table_tsv = self.taxon_table_tsv_dict[self.beta_diversity_feature_level]

        rooted_tree_qza = Phylogeny(self.settings).main(
            seq_qza=self.feature_sequence_qza)  # no need to use taxonomy-labeled feature sequences

        QiimeBetaDiversity(self.settings).main(
            feature_table_tsv=feature_table_tsv,
            rooted_tree_qza=rooted_tree_qza,
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
        table_tsv_dict.update(self.picrust2_table_tsv_dict)
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
            colors=self.colors,
            p_value=self.differential_abundance_p_value)

    def collect_log_files(self):
        makedirs(f'{self.outdir}/log', exist_ok=True)
        cmd = f'mv "{self.outdir}"/*.log "{self.outdir}"/log/'
        self.call(cmd)
