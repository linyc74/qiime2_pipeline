from .template import Processor, Settings


class CutadaptTrimPaired(Processor):

    untrimmed_reads_qza: str

    trimmed_reads_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            untrimmed_reads_qza: str) -> str:

        self.untrimmed_reads_qza = untrimmed_reads_qza

        self.set_trimmed_reads_qza()
        self.execute()

        return self.trimmed_reads_qza

    def set_trimmed_reads_qza(self):
        self.trimmed_reads_qza = f'{self.workdir}/cutadapt-trimmed-reads.qza'

    def execute(self):
        lines = [
            'qiime cutadapt trim-paired',
            f'--i-demultiplexed-sequences {self.untrimmed_reads_qza}',
            f'--p-cores {self.threads}',
            f'--o-trimmed-sequences {self.trimmed_reads_qza}',
        ]
        self.call(' \\\n  '.join(lines))
