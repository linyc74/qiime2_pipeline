import os
import gzip
from typing import IO, List, Tuple
from .tools import get_files, rev_comp
from .template import Processor, Settings


class Concat(Processor):

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
        self.writer = open(self.fq_out, 'wt')

    def concat(self):
        i = 0
        for line1, line2 in zip(self.reader1, self.reader2):
            i += 1

            is_seq = i % 4 == 2
            is_qual = i % 4 == 0
            l1 = line1.strip()
            l2 = line2.strip()

            if is_seq:
                new = f'{l1}{rev_comp(l2)}\n'
            elif is_qual:
                new = f'{l1}{l2[::-1]}\n'
            else:
                new = f'{l1}\n'

            self.writer.write(new)

    def close_files(self):
        self.reader1.close()
        self.reader2.close()
        self.writer.close()


class BatchConcat(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str

    sample_names: List[str]
    out_fq_dir: str
    out_fq_suffix: str = '.fq'

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.concat = Concat(self.settings).main

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str) -> Tuple[str, str]:

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix

        self.set_sample_names()
        self.set_out_fq_dir()
        for name in self.sample_names:
            self.process_one_pair(name)

        return self.out_fq_dir, self.out_fq_suffix

    def set_sample_names(self):
        files = get_files(
            source=self.fq_dir,
            endswith=self.fq1_suffix,
            isfullpath=False)
        n_char = len(self.fq1_suffix)
        self.sample_names = [f[:-n_char] for f in files]

    def set_out_fq_dir(self):
        self.out_fq_dir = f'{self.workdir}/concat_fastqs'
        os.makedirs(self.out_fq_dir, exist_ok=True)

    def process_one_pair(self, name: str):
        fq1 = f'{self.fq_dir}/{name}{self.fq1_suffix}'
        fq2 = f'{self.fq_dir}/{name}{self.fq2_suffix}'

        fq = self.concat(fq1, fq2)

        os.rename(
            fq, f'{self.out_fq_dir}/{name}{self.out_fq_suffix}'
        )
