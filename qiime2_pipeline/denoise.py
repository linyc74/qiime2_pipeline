from typing import Tuple
from .template import Processor, Settings


class Dada2Paired(Processor):

    TRIM_LEFT = 0
    TRUNCATE_LENGTH = 0
    MIN_OVERLAP = 4

    demultiplexed_seq_qza: str
    truncate_length: int

    feature_sequence_qza: str
    feature_sequence_fa: str
    feature_table_qza: str
    feature_table_tsv: str
    denoising_stats_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            demultiplexed_seq_qza: str) -> Tuple[str, str]:

        self.demultiplexed_seq_qza = demultiplexed_seq_qza

        self.set_output_paths()
        self.execute()
        self.export()

        return self.feature_sequence_qza, self.feature_table_qza

    def set_output_paths(self):
        self.feature_sequence_qza = f'{self.workdir}/dada2-feature-sequence.qza'
        self.feature_sequence_fa = f'{self.outdir}/dada2-feature-sequence.fa'
        self.feature_table_qza = f'{self.workdir}/dada2-feature-table.qza'
        self.feature_table_tsv = f'{self.outdir}/dada2-feature-table.tsv'
        self.denoising_stats_qza = f'{self.workdir}/dada2-stats.qza'

    def execute(self):
        lines = [
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
        ]
        self.call(' \\\n  '.join(lines))

    def export(self):
        ExportFeatureSequence(self.settings).main(
            feature_sequence_qza=self.feature_sequence_qza,
            output_fa=self.feature_sequence_fa
        )

        ExportFeatureTable(self.settings).main(
            feature_table_qza=self.feature_table_qza,
            output_tsv=self.feature_table_tsv
        )


class ExportFeatureTable(Processor):

    feature_table_qza: str
    output_tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_table_qza: str, output_tsv: str):
        self.feature_table_qza = feature_table_qza
        self.output_tsv = output_tsv

        self.qza_to_biom()
        self.biom_to_tsv()

    def qza_to_biom(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.feature_table_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

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

    def main(
            self,
            feature_sequence_qza: str,
            output_fa: str):

        self.feature_sequence_qza = feature_sequence_qza
        self.output_fa = output_fa

        self.qza_to_fa()
        self.move_fa()

    def qza_to_fa(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.feature_sequence_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def move_fa(self):
        self.call(f'mv {self.workdir}/dna-sequences.fasta {self.output_fa}')


class VsearchPaired(Processor):

    demultiplexed_seq_qza: str
    truncate_length: int

    joined_sequences: str
    representative_seq_qza: str
    table_qza: str
    stats_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            demultiplexed_seq_qza: str,
            truncate_length: int) -> Tuple[str, str]:

        self.demultiplexed_seq_qza = demultiplexed_seq_qza
        self.truncate_length = truncate_length

        self.set_output_paths()
        self.execute()

        return self.representative_seq_qza, self.table_qza

    def set_output_paths(self):
        self.joined_sequences = f'{self.workdir}/vsearch-joined-sequences.qza'
        self.representative_seq_qza = f'{self.workdir}/vsearch-representative-sequences.qza'
        self.table_qza = f'{self.workdir}/vsearch-table.qza'
        self.stats_qza = f'{self.workdir}/vsearch-stats.qza'

    def execute(self):
        lines = [
            'qiime vsearch join-pairs',
            f'--i-demultiplexed-seqs {self.demultiplexed_seq_qza}',
            f'--o-joined-sequences {self.joined_sequences}',
        ]
        self.call(' \\\n  '.join(lines))

        lines = [
            ''
        ]
