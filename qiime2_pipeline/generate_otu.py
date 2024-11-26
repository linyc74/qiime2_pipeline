from typing import Tuple
from .template import Processor
from .importing import ImportSingleEndFastq


class GenerateOTU(Processor):

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

        self.clustered_table_qza = f'{self.workdir}/otu-feature-table.qza'
        self.clustered_sequence_qza = f'{self.workdir}/otu-feature-sequence.qza'
        log = f'{self.outdir}/qiime-vsearch-cluster-features-de-novo.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime vsearch cluster-features-de-novo',
            f'--i-sequences {self.feature_sequence_qza}',
            f'--i-table {self.feature_table_qza}',
            f'--p-perc-identity {self.identity}',
            f'--p-threads {self.threads}',
            f'--o-clustered-table {self.clustered_table_qza}',
            f'--o-clustered-sequences {self.clustered_sequence_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)

        return self.clustered_table_qza, self.clustered_sequence_qza


class GenerateNanoporeOTU(Processor):

    sample_sheet: str
    fq_dir: str
    fq_suffix: str
    identity: float

    single_end_seq_qza: str

    dereplicated_table_qza: str
    dereplicated_sequence_qza: str

    clustered_table_qza: str
    clustered_sequence_qza: str

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq_suffix: str,
            identity: float) -> Tuple[str, str]:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq_suffix = fq_suffix
        self.identity = identity

        self.single_end_seq_qza = ImportSingleEndFastq(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq_suffix=self.fq_suffix)

        self.dereplicate_sequences()

        self.clustered_table_qza, self.clustered_sequence_qza = GenerateOTU(self.settings).main(
            feature_table_qza=self.dereplicated_table_qza,
            feature_sequence_qza=self.dereplicated_sequence_qza,
            identity=self.identity)

        return self.clustered_table_qza, self.clustered_sequence_qza

    def dereplicate_sequences(self):
        self.dereplicated_table_qza = f'{self.workdir}/dereplicated-table.qza'
        self.dereplicated_sequence_qza = f'{self.workdir}/dereplicated-sequence.qza'
        log = f'{self.outdir}/qiime-vsearch-dereplicate-sequences.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime vsearch dereplicate-sequences',
            f'--i-sequences {self.single_end_seq_qza}',
            f'--o-dereplicated-table {self.dereplicated_table_qza}',
            f'--o-dereplicated-sequences {self.dereplicated_sequence_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)
