from typing import Tuple
from .template import Processor


class Dada2(Processor):

    TRIM_LEFT = 0  # number of 5' bases to be clipped
    TRUNCATE_LENGTH = 0  # min read length

    demultiplexed_seq_qza: str
    max_expected_error_bases: float

    feature_sequence_qza: str
    feature_table_qza: str
    denoising_stats_qza: str

    def main(
            self,
            demultiplexed_seq_qza: str,
            max_expected_error_bases: float) -> Tuple[str, str]:

        self.demultiplexed_seq_qza = demultiplexed_seq_qza
        self.max_expected_error_bases = max_expected_error_bases

        self.set_output_paths()
        self.execute()
        self.export_stats()

        return self.feature_table_qza, self.feature_sequence_qza

    def set_output_paths(self):
        self.feature_sequence_qza = f'{self.workdir}/dada2-feature-sequence.qza'
        self.feature_table_qza = f'{self.workdir}/dada2-feature-table.qza'
        self.denoising_stats_qza = f'{self.workdir}/dada2-stats.qza'

    def execute(self):
        pass

    def export_stats(self):
        out = f'{self.workdir}/dada2'
        log = f'{self.outdir}/qiime-tools-export.log'

        cmd = self.CMD_LINEBREAK.join([
            'qiime tools export',
            f'--input-path {self.denoising_stats_qza}',
            f'--output-path {out}',
            f'1>> {log}',
            f'2>> {log}'
        ])
        self.call(cmd)

        self.call(f'mv {out}/stats.tsv {self.outdir}/dada2-stats.tsv')


class Dada2SingleEnd(Dada2):

    def execute(self):
        log = f'{self.outdir}/qiime-dada2-denoise-single.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime dada2 denoise-single',
            f'--i-demultiplexed-seqs {self.demultiplexed_seq_qza}',
            f'--p-trim-left {self.TRIM_LEFT}',
            f'--p-trunc-len {self.TRUNCATE_LENGTH}',
            f'--p-max-ee {self.max_expected_error_bases}',
            f'--p-n-threads {self.threads}',
            f'--o-representative-sequences {self.feature_sequence_qza}',
            f'--o-table {self.feature_table_qza}',
            f'--o-denoising-stats {self.denoising_stats_qza}',
            f'1>> {log}',
            f'2>> {log}'
        ])
        self.call(cmd)


class Dada2PairedEnd(Dada2):

    MIN_OVERLAP = 12

    def execute(self):
        log = f'{self.outdir}/qiime-dada2-denoise-paired.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime dada2 denoise-paired',
            f'--i-demultiplexed-seqs {self.demultiplexed_seq_qza}',
            f'--p-trim-left-f {self.TRIM_LEFT}',
            f'--p-trim-left-r {self.TRIM_LEFT}',
            f'--p-trunc-len-f {self.TRUNCATE_LENGTH}',
            f'--p-trunc-len-r {self.TRUNCATE_LENGTH}',
            f'--p-max-ee-f {self.max_expected_error_bases}',
            f'--p-max-ee-r {self.max_expected_error_bases}',
            f'--p-min-overlap {self.MIN_OVERLAP}',
            f'--p-n-threads {self.threads}',
            f'--o-representative-sequences {self.feature_sequence_qza}',
            f'--o-table {self.feature_table_qza}',
            f'--o-denoising-stats {self.denoising_stats_qza}',
            f'1>> {log}',
            f'2>> {log}'
        ])
        self.call(cmd)
