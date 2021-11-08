from typing import Tuple
from .template import Processor, Settings


class Dada2Paired(Processor):

    TRIM_LEFT = 0

    demultiplexed_seq_qza: str
    truncate_length: int

    representative_seq_qza: str
    table_qza: str
    denoising_stats_qza: str

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
        self.representative_seq_qza = f'{self.workdir}/dada2-representative-sequences.qza'
        self.table_qza = f'{self.workdir}/dada2-table.qza'
        self.denoising_stats_qza = f'{self.workdir}/dada2-stats.qza'

    def execute(self):
        lines = [
            'qiime dada2 denoise-paired',
            f'--i-demultiplexed-seqs {self.demultiplexed_seq_qza}',
            f'--p-trim-left-f {self.TRIM_LEFT}',
            f'--p-trim-left-r {self.TRIM_LEFT}',
            f'--p-trunc-len-f {self.truncate_length}',
            f'--p-trunc-len-r {self.truncate_length}',
            f'--p-n-threads {self.threads}',
            f'--o-representative-sequences {self.representative_seq_qza}',
            f'--o-table {self.table_qza}',
            f'--o-denoising-stats {self.denoising_stats_qza}',
        ]
        self.call(' \\\n  '.join(lines))
