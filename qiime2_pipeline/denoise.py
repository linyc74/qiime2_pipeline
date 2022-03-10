from typing import Tuple
from .template import Processor


class Dada2(Processor):

    TRIM_LEFT = 0
    TRUNCATE_LENGTH = 0

    demultiplexed_seq_qza: str

    feature_sequence_qza: str
    feature_table_qza: str
    denoising_stats_qza: str

    def set_output_paths(self):
        self.feature_sequence_qza = f'{self.workdir}/dada2-feature-sequence.qza'
        self.feature_table_qza = f'{self.workdir}/dada2-feature-table.qza'
        self.denoising_stats_qza = f'{self.workdir}/dada2-stats.qza'


class Dada2SingleEnd(Dada2):

    def main(
            self,
            demultiplexed_seq_qza: str) -> Tuple[str, str]:

        self.demultiplexed_seq_qza = demultiplexed_seq_qza

        self.set_output_paths()
        self.execute()

        return self.feature_table_qza, self.feature_sequence_qza

    def execute(self):
        cmd = self.CMD_LINEBREAK.join([
            'qiime dada2 denoise-single',
            f'--i-demultiplexed-seqs {self.demultiplexed_seq_qza}',
            f'--p-trim-left {self.TRIM_LEFT}',
            f'--p-trunc-len {self.TRUNCATE_LENGTH}',
            f'--p-n-threads {self.threads}',
            f'--o-representative-sequences {self.feature_sequence_qza}',
            f'--o-table {self.feature_table_qza}',
            f'--o-denoising-stats {self.denoising_stats_qza}',
        ])
        self.call(cmd)


class Dada2PairedEnd(Dada2):

    MIN_OVERLAP = 12

    def main(
            self,
            demultiplexed_seq_qza: str) -> Tuple[str, str]:

        self.demultiplexed_seq_qza = demultiplexed_seq_qza

        self.set_output_paths()
        self.execute()

        return self.feature_table_qza, self.feature_sequence_qza

    def execute(self):
        cmd = self.CMD_LINEBREAK.join([
            'qiime dada2 denoise-paired',
            f'--i-demultiplexed-seqs {self.demultiplexed_seq_qza}',
            f'--p-trim-left-f {self.TRIM_LEFT}',
            f'--p-trim-left-r {self.TRIM_LEFT}',
            f'--p-trunc-len-f {self.TRUNCATE_LENGTH}',
            f'--p-trunc-len-r {self.TRUNCATE_LENGTH}',
            f'--p-min-overlap {self.MIN_OVERLAP}',
            f'--p-n-threads {self.threads}',
            f'--o-representative-sequences {self.feature_sequence_qza}',
            f'--o-table {self.feature_table_qza}',
            f'--o-denoising-stats {self.denoising_stats_qza}',
        ])
        self.call(cmd)
