import gzip
import pandas as pd
from os.path import basename
from typing import List, Dict, Union, Optional
from .tools import get_files
from .template import Processor


class RawReadCounts(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: Optional[str]

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: Optional[str]):

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix

        if self.fq2_suffix is None:
            SingleEnd(self.settings).main(
                fq_dir=self.fq_dir,
                fq_suffix=self.fq1_suffix)
        else:
            PairedEnd(self.settings).main(
                fq_dir=self.fq_dir,
                fq1_suffix=self.fq1_suffix,
                fq2_suffix=self.fq2_suffix)


class PairedEnd(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str

    fq1s: List[str]
    fq2s: List[str]
    data: List[Dict[str, Union[str, int]]]

    def main(self, fq_dir: str, fq1_suffix: str, fq2_suffix: str):
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix

        self.set_fqs()
        self.read_fqs()
        self.save_csv()

    def set_fqs(self):
        self.fq1s = get_files(
            source=self.fq_dir,
            endswith=self.fq1_suffix,
            isfullpath=True)
        self.fq2s = get_files(
            source=self.fq_dir,
            endswith=self.fq1_suffix,
            isfullpath=True)
        assert len(self.fq1s) == len(self.fq2s), f'The number of read 1 fastq files is not equal to the number of read 2 fastq files'

    def read_fqs(self):
        self.data = []
        for fq1, fq2 in zip(self.fq1s, self.fq2s):
            self.data.append({
                'Sample ID': basename(fq1)[:-len(self.fq1_suffix)],
                'Count (R1)': count_reads(fq1),
                'Count (R2)': count_reads(fq2),
            })

    def save_csv(self):
        pd.DataFrame(self.data).to_csv(f'{self.outdir}/raw-read-counts.csv', index=False)


class SingleEnd(Processor):

    fq_dir: str
    fq_suffix: str

    fqs: List[str]
    data: List[Dict[str, Union[str, int]]]

    def main(self, fq_dir: str, fq_suffix: str):
        self.fq_dir = fq_dir
        self.fq_suffix = fq_suffix

        self.set_fqs()
        self.read_fqs()
        self.save_csv()

    def set_fqs(self):
        self.fqs = get_files(
            source=self.fq_dir,
            endswith=self.fq_suffix,
            isfullpath=True)

    def read_fqs(self):
        self.data = []
        for fq in self.fqs:
            self.data.append({
                'Sample ID': basename(fq)[:-len(self.fq_suffix)],
                'Count': count_reads(fq)
            })

    def save_csv(self):
        pd.DataFrame(self.data).to_csv(f'{self.outdir}/raw-read-counts.csv', index=False)


def count_reads(fq: str) -> int:
    fh = gzip.open(fq) if fq.endswith('.gz') else open(fq)

    i = 0
    for _ in fh:
        i += 1

    fh.close()

    assert i % 4 == 0
    return i // 4
