import os
from typing import List, Union
from .template import Settings
from .qiime2_pipeline import Qiime2Pipeline


class Main:

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str
    nb_classifier_qza: str
    paired_end_mode: str
    group_keywords: List[str]
    otu_identity: float
    skip_otu: bool
    classifier_reads_per_batch: Union[int, str]

    settings: Settings

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            nb_classifier_qza: str,
            paired_end_mode: str,
            group_keywords: str,
            otu_identity: float,
            skip_otu: bool,
            classifier_reads_per_batch: int,
            outdir: str,
            threads: int,
            debug: bool):

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.nb_classifier_qza = nb_classifier_qza
        self.paired_end_mode = paired_end_mode
        self.group_keywords = [] if group_keywords == 'None' \
            else group_keywords.split(',')
        self.otu_identity = float(otu_identity)
        self.skip_otu = skip_otu
        self.classifier_reads_per_batch = classifier_reads_per_batch

        self.settings = Settings(
            workdir='./qiime2_pipeline_workdir',
            outdir=outdir,
            threads=int(threads),
            debug=debug,
            mock=False)

        for d in [self.settings.workdir, self.settings.outdir]:
            os.makedirs(d, exist_ok=True)

        Qiime2Pipeline(self.settings).main(
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            nb_classifier_qza=self.nb_classifier_qza,
            paired_end_mode=self.paired_end_mode,
            group_keywords=self.group_keywords,
            otu_identity=self.otu_identity,
            skip_otu=self.skip_otu,
            classifier_reads_per_batch=self.classifier_reads_per_batch)
