import os
import shutil
from typing import List, Optional
from .template import Settings
from .tools import get_temp_path
from .qiime2_pipeline import Qiime2Pipeline


class Main:

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
    invert_colors: bool

    settings: Settings

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            nb_classifier_qza: str,
            paired_end_mode: str,
            otu_identity: float,
            skip_otu: bool,
            classifier_reads_per_batch: int,
            alpha_metrics: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int,
            max_expected_error_bases: float,
            heatmap_read_fraction: float,
            n_taxa_barplot: int,
            colormap: str,
            invert_colors: bool,
            outdir: str,
            threads: int,
            debug: bool):

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = None if fq2_suffix == 'None' else fq2_suffix
        self.nb_classifier_qza = nb_classifier_qza
        self.paired_end_mode = paired_end_mode
        self.otu_identity = float(otu_identity)
        self.skip_otu = skip_otu
        self.classifier_reads_per_batch = classifier_reads_per_batch
        self.alpha_metrics = [] if alpha_metrics == 'all' else alpha_metrics.split(',')
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime
        self.max_expected_error_bases = max_expected_error_bases
        self.heatmap_read_fraction = heatmap_read_fraction
        self.n_taxa_barplot = n_taxa_barplot
        self.colormap = colormap
        self.invert_colors = invert_colors

        self.settings = Settings(
            workdir=get_temp_path(prefix='./qiime2_pipeline_workdir_'),
            outdir=outdir,
            threads=int(threads),
            debug=debug,
            mock=False)

        self.makedirs()
        self.run_pipeline()
        self.remove_workdir()

    def makedirs(self):
        for d in [self.settings.workdir, self.settings.outdir]:
            os.makedirs(d, exist_ok=True)

    def run_pipeline(self):
        Qiime2Pipeline(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            nb_classifier_qza=self.nb_classifier_qza,
            paired_end_mode=self.paired_end_mode,
            otu_identity=self.otu_identity,
            skip_otu=self.skip_otu,
            classifier_reads_per_batch=self.classifier_reads_per_batch,
            alpha_metrics=self.alpha_metrics,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime,
            max_expected_error_bases=self.max_expected_error_bases,
            heatmap_read_fraction=self.heatmap_read_fraction,
            n_taxa_barplot=self.n_taxa_barplot,
            colormap=self.colormap,
            invert_colors=self.invert_colors)

    def remove_workdir(self):
        if not self.settings.debug:
            shutil.rmtree(self.settings.workdir)
