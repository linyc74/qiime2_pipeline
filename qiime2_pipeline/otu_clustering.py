from typing import Tuple
from .template import Processor


class Vsearch(Processor):

    feature_table_qza: str
    feature_sequence_qza: str
    identity: float

    clustered_table_qza: str
    clustered_sequence_qza: str

    def main(
            self,
            feature_table_qza: str,
            feature_sequence_qza: str,
            identity: float) -> Tuple[str, str]:

        self.feature_table_qza = feature_table_qza
        self.feature_sequence_qza = feature_sequence_qza
        self.identity = identity

        self.set_output_paths()
        self.execute()

        return self.clustered_table_qza, self.clustered_sequence_qza

    def set_output_paths(self):
        self.clustered_table_qza = f'{self.workdir}/vsearch-feature-table.qza'
        self.clustered_sequence_qza = f'{self.workdir}/vsearch-feature-sequence.qza'

    def execute(self):
        cmd = self.CMD_LINEBREAK.join([
            'qiime vsearch cluster-features-de-novo',
            f'--i-sequences {self.feature_sequence_qza}',
            f'--i-table {self.feature_table_qza}',
            f'--p-perc-identity {self.identity}',
            f'--p-threads {self.threads}',
            f'--o-clustered-table {self.clustered_table_qza}',
            f'--o-clustered-sequences {self.clustered_sequence_qza}',
        ])
        self.call(cmd)
