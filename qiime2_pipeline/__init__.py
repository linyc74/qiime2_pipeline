import os
import shutil
from .template import Settings
from .utils import get_temp_path
from .qiime2_pipeline import Qiime2Pipeline


def main(
        self,
        sample_sheet: str,
        fq_dir: str,
        fq1_suffix: str,
        fq2_suffix: str,
        pacbio: bool,
        paired_end_mode: str,
        otu_identity: float,
        skip_otu: bool,
        feature_classifier: str,
        nb_classifier_qza: str,
        classifier_reads_per_batch: int,
        reference_sequence_qza: str,
        reference_taxonomy_qza: str,
        alpha_metrics: str,
        clip_r1_5_prime: int,
        clip_r2_5_prime: int,
        max_expected_error_bases: float,
        heatmap_read_fraction: float,
        n_taxa_barplot: int,
        beta_diversity_feature_level: str,
        colormap: str,
        invert_colors: bool,
        publication_figure: bool,
        skip_differential_abundance: bool,
        run_picrust2: bool,
        outdir: str,
        threads: int,
        debug: bool):

    settings = Settings(
        workdir=get_temp_path(prefix='./qiime2_pipeline_workdir_'),
        outdir=outdir,
        threads=int(threads),
        debug=debug,
        mock=False,
        for_publication=publication_figure)

    for d in [settings.workdir, outdir]:
        os.makedirs(d, exist_ok=True)

    Qiime2Pipeline(settings).main(
        sample_sheet=sample_sheet,
        fq_dir=fq_dir,
        fq1_suffix=fq1_suffix,
        fq2_suffix=None if fq2_suffix.lower() == 'none' else fq2_suffix,
        pacbio=pacbio,
        paired_end_mode=paired_end_mode,
        otu_identity=otu_identity,
        skip_otu=skip_otu,
        feature_classifier=feature_classifier,
        nb_classifier_qza=None if nb_classifier_qza.lower() == 'none' else nb_classifier_qza,
        classifier_reads_per_batch=classifier_reads_per_batch,
        reference_sequence_qza=None if reference_sequence_qza.lower() == 'none' else reference_sequence_qza,
        reference_taxonomy_qza=None if reference_taxonomy_qza.lower() == 'none' else reference_taxonomy_qza,
        alpha_metrics=[] if alpha_metrics == 'all' else alpha_metrics.split(','),
        clip_r1_5_prime=clip_r1_5_prime,
        clip_r2_5_prime=clip_r2_5_prime,
        max_expected_error_bases=max_expected_error_bases,
        heatmap_read_fraction=heatmap_read_fraction,
        n_taxa_barplot=n_taxa_barplot,
        beta_diversity_feature_level=beta_diversity_feature_level,
        colormap=colormap,
        invert_colors=invert_colors,
        skip_differential_abundance=skip_differential_abundance,
        run_picrust2=run_picrust2)

    if not debug:
        shutil.rmtree(settings.workdir)
