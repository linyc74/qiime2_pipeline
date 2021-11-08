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
            nb_classifier_qza: str) -> str:

        self.representative_seq_qza = representative_seq_qza
        self.nb_classifier_qza = nb_classifier_qza

        self.set_taxonomy_qza()
        self.execute()

        return self.taxonomy_qza

    def set_taxonomy_qza(self):
        self.taxonomy_qza = f'{self.workdir}/taxonomy.qza'

    def execute(self):
        lines = [
            'qiime feature-classifier classify-sklearn',
            f'--i-classifier {self.nb_classifier_qza}',
            f'--i-reads {self.representative_seq_qza}',
            f'--p-n-jobs {self.threads}',
            f'--o-classification {self.taxonomy_qza}',
        ]
        self.call(' \\\n  '.join(lines))
