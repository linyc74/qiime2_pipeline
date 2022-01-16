import gzip
from .template import Processor, Settings
from typing import IO


class ConcatRead1Read2(Processor):

    fq1: str
    fq2: str

    reader1: IO
    reader2: IO
    writer: IO
    fq_out: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            fq1: str,
            fq2: str) -> str:

        self.fq1 = fq1
        self.fq2 = fq2

        self.open_files()
        self.concat()
        self.close_files()

        return self.fq_out

    def open_files(self):
        self.reader1 = gzip.open(self.fq1, 'rt')
        self.reader2 = gzip.open(self.fq2, 'rt')
        self.fq_out = f'{self.workdir}/concat.fq'
        self.writer = open(self.fq_out, 'w')

    def concat(self):
        i = 0
        for line1, line2 in zip(self.reader1, self.reader2):
            i += 1

            is_seq = i % 4 == 2
            is_qual = i % 4 == 0

            if is_seq or is_qual:
                self.writer.write(f'{line1.strip()}{line2}')
            else:
                self.writer.write(line1)

    def close_files(self):
        self.reader1.close()
        self.reader2.close()
        self.writer.close()
