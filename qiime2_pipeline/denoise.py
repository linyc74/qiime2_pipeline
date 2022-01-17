from typing import Tuple
from os.path import basename
from .template import Processor, Settings


class Dada2(Processor):

    TRIM_LEFT = 0
    TRUNCATE_LENGTH = 0

    demultiplexed_seq_qza: str

    feature_sequence_qza: str
    feature_table_qza: str
    denoising_stats_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def set_output_paths(self):
        self.feature_sequence_qza = f'{self.workdir}/dada2-feature-sequence.qza'
        self.feature_table_qza = f'{self.workdir}/dada2-feature-table.qza'
        self.denoising_stats_qza = f'{self.workdir}/dada2-stats.qza'

    def export(self):
        ExportFeatureSequence(self.settings).main(
            feature_sequence_qza=self.feature_sequence_qza)

        ExportFeatureTable(self.settings).main(
            feature_table_qza=self.feature_table_qza)


class Dada2SingleEnd(Dada2):

    def main(
            self,
            demultiplexed_seq_qza: str) -> Tuple[str, str]:

        self.demultiplexed_seq_qza = demultiplexed_seq_qza

        self.set_output_paths()
        self.execute()
        self.export()

        return self.feature_sequence_qza, self.feature_table_qza

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
        self.export()

        return self.feature_sequence_qza, self.feature_table_qza

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


class ExportFeatureTable(Processor):

    feature_table_qza: str
    output_tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_table_qza: str):
        self.feature_table_qza = feature_table_qza

        self.qza_to_biom()
        self.set_output_tsv()
        self.biom_to_tsv()

    def qza_to_biom(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.feature_table_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def set_output_tsv(self):
        fname = basename(self.feature_table_qza)[:-len('.qza')] + '.tsv'
        self.output_tsv = f'{self.outdir}/{fname}'

    def biom_to_tsv(self):
        lines = [
            'biom convert --to-tsv',
            f'-i {self.workdir}/feature-table.biom',
            f'-o {self.output_tsv}'
        ]
        self.call(' \\\n  '.join(lines))


class ExportFeatureSequence(Processor):

    feature_sequence_qza: str
    output_fa: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_sequence_qza: str):

        self.feature_sequence_qza = feature_sequence_qza

        self.qza_to_fa()
        self.set_output_fa()
        self.move_fa()

    def qza_to_fa(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.feature_sequence_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def set_output_fa(self):
        fname = basename(self.feature_sequence_qza)[:-len('.qza')] + '.fa'
        self.output_fa = f'{self.outdir}/{fname}'

    def move_fa(self):
        self.call(f'mv {self.workdir}/dna-sequences.fasta {self.output_fa}')
