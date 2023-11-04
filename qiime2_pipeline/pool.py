import os
import pandas as pd
from typing import List, Tuple
from .template import Processor, Settings


class Pool(Processor):

    fq1: str
    fq2: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            fq1: str,
            fq2: str) -> str:

        self.fq1 = fq1
        self.fq2 = fq2

        f = f'{self.workdir}/pool.fq.gz'
        self.call(f'cat {self.fq1} {self.fq2} > {f}')
        self.call(f'gunzip {f}')

        return f[:-len('.gz')]


class BatchPool(Processor):

    OUT_FQ_DIRNAME = 'pool_fastqs'
    OUT_FQ_SUFFIX = '.fq'

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str

    sample_names: List[str]
    out_fq_dir: str

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str) -> Tuple[str, str]:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix

        self.set_sample_names()
        self.set_out_fq_dir()
        for name in self.sample_names:
            self.process_one_pair(name)

        return self.out_fq_dir, self.OUT_FQ_SUFFIX

    def set_sample_names(self):
        self.sample_names = pd.read_csv(self.sample_sheet, index_col=0).index.tolist()

    def set_out_fq_dir(self):
        self.out_fq_dir = f'{self.workdir}/{self.OUT_FQ_DIRNAME}'
        os.makedirs(self.out_fq_dir, exist_ok=True)

    def process_one_pair(self, name: str):
        fq1 = f'{self.fq_dir}/{name}{self.fq1_suffix}'
        fq2 = f'{self.fq_dir}/{name}{self.fq2_suffix}'

        fq = Pool(self.settings).main(fq1, fq2)

        os.rename(
            fq, f'{self.out_fq_dir}/{name}{self.OUT_FQ_SUFFIX}'
        )
