import gzip
import pandas as pd
from os.path import basename
from typing import List, Dict, Union, Optional
from .template import Processor


class RawReadCounts(Processor):

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: Optional[str]

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: Optional[str]):

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix

        if self.fq2_suffix is None:
            SingleEnd(self.settings).main(
                sample_sheet=self.sample_sheet,
                fq_dir=self.fq_dir,
                fq_suffix=self.fq1_suffix)
        else:
            PairedEnd(self.settings).main(
                sample_sheet=self.sample_sheet,
                fq_dir=self.fq_dir,
                fq1_suffix=self.fq1_suffix,
                fq2_suffix=self.fq2_suffix)


class PairedEnd(Processor):

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str

    fq1s: List[str]
    fq2s: List[str]
    data: List[Dict[str, Union[str, int]]]

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str):

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix

        self.set_fqs()
        self.read_fqs()
        self.save_csv()

    def set_fqs(self):
        df = pd.read_csv(self.sample_sheet, index_col=0)
        self.fq1s = [
            f'{self.fq_dir}/{sample_id}{self.fq1_suffix}'
            for sample_id in df.index
        ]
        self.fq2s = [
            f'{self.fq_dir}/{sample_id}{self.fq2_suffix}'
            for sample_id in df.index
        ]

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

    sample_sheet: str
    fq_dir: str
    fq_suffix: str

    fqs: List[str]
    data: List[Dict[str, Union[str, int]]]

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq_suffix: str):

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq_suffix = fq_suffix

        self.set_fqs()
        self.read_fqs()
        self.save_csv()

    def set_fqs(self):
        df = pd.read_csv(self.sample_sheet, index_col=0)
        self.fqs = [
            f'{self.fq_dir}/{sample_id}{self.fq_suffix}'
            for sample_id in df.index
        ]

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
