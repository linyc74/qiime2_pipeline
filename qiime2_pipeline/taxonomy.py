from os.path import basename
from .template import Processor, Settings


class FeatureClassifier(Processor):

    representative_seq_qza: str
    nb_classifier_qza: str

    taxonomy_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            representative_seq_qza: str,
            nb_classifier_qza: str):

        self.representative_seq_qza = representative_seq_qza
        self.nb_classifier_qza = nb_classifier_qza

        self.execute()
        self.export()

    def execute(self):
        self.taxonomy_qza = f'{self.workdir}/taxonomy.qza'
        lines = [
            'qiime feature-classifier classify-sklearn',
            f'--i-classifier {self.nb_classifier_qza}',
            f'--i-reads {self.representative_seq_qza}',
            f'--p-n-jobs {self.threads}',
            f'--o-classification {self.taxonomy_qza}',
        ]
        self.call(' \\\n  '.join(lines))

    def export(self):
        ExportTaxonomy(self.settings).main(taxonomy_qza=self.taxonomy_qza)


class ExportTaxonomy(Processor):

    taxonomy_qza: str
    output_tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, taxonomy_qza: str):

        self.taxonomy_qza = taxonomy_qza

        self.set_output_tsv()
        self.qza_to_tsv()
        self.move_tsv()

    def qza_to_tsv(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.taxonomy_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def set_output_tsv(self):
        fname = basename(self.taxonomy_qza)[:-len('.qza')] + '.tsv'
        self.output_tsv = f'{self.outdir}/{fname}'

    def move_tsv(self):
        self.call(f'mv {self.workdir}/taxonomy.tsv {self.output_tsv}')
