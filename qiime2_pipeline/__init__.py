import os
from .template import Settings
from .qiime2_pipeline import Qiime2Pipeline


class Main:

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str
    nb_classifier_qza: str

    settings: Settings

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            nb_classifier_qza: str,
            outdir: str,
            threads: str,
            debug: bool):

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.nb_classifier_qza = nb_classifier_qza

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
            nb_classifier_qza=self.nb_classifier_qza)
