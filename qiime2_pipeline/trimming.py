import os
import pandas as pd
from os.path import basename
from typing import Tuple, List


from .template import Processor, Settings


class TrimGalorePairedEnd(Processor):

    QUALITY = 20
    LENGTH = 20
    MAX_N = 0
    CUTADAPT_TOTAL_CORES = 2
    # According to the help message of trim_galore, 2 cores for cutadapt -> actually up to 9 cores

    fq1: str
    fq2: str
    clip_r1_5_prime: int
    clip_r2_5_prime: int

    out_fq1: str
    out_fq2: str

    def main(
            self,
            fq1: str,
            fq2: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int) -> Tuple[str, str]:

        self.fq1 = fq1
        self.fq2 = fq2
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        self.execute()
        self.move_fastqc_report()
        self.set_out_fq1()
        self.set_out_fq2()

        return self.out_fq1, self.out_fq2

    def execute(self):
        args = [
            'trim_galore',
            '--paired',
            f'--quality {self.QUALITY}',
            '--phred33',
            f'--cores {self.CUTADAPT_TOTAL_CORES}',
            f'--fastqc_args "--threads {self.threads}"',
            '--illumina',
            f'--length {self.LENGTH}',
            f'--max_n {self.MAX_N}',
            '--trim-n',
            '--gzip',
            f'--output_dir {self.workdir}'
        ]

        if self.clip_r1_5_prime > 0:
            args.append(f'--clip_R1 {self.clip_r1_5_prime}')

        if self.clip_r2_5_prime > 0:
            args.append(f'--clip_R2 {self.clip_r2_5_prime}')

        log = f'{self.outdir}/trim_galore.log'
        args += [
            self.fq1,
            self.fq2,
            f'1>> {log} 2>> {log}'
        ]

        self.call(self.CMD_LINEBREAK.join(args))

    def move_fastqc_report(self):
        dstdir = f'{self.outdir}/fastqc'
        os.makedirs(dstdir, exist_ok=True)
        for suffix in [
            'fastqc.html',
            'fastqc.zip',
            'trimming_report.txt'
        ]:
            self.call(f'mv {self.workdir}/*{suffix} {dstdir}/')

    def set_out_fq1(self):
        f = basename(self.fq1)
        f = self.__strip_file_extension(f)
        self.out_fq1 = f'{self.workdir}/{f}_val_1.fq.gz'

    def set_out_fq2(self):
        f = basename(self.fq2)
        f = self.__strip_file_extension(f)
        self.out_fq2 = f'{self.workdir}/{f}_val_2.fq.gz'

    def __strip_file_extension(self, f):
        for suffix in [
            '.fq',
            '.fq.gz',
            '.fastq',
            '.fastq.gz',
        ]:
            if f.endswith(suffix):
                f = f[:-len(suffix)]  # strip suffix
        return f


class BatchTrimGalorePairedEnd(Processor):

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str
    clip_r1_5_prime: int
    clip_r2_5_prime: int

    sample_names: List[str]
    out_fq_dir: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.trim_galore = TrimGalorePairedEnd(self.settings).main

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int) -> str:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        self.set_sample_names()
        self.set_out_fq_dir()
        for name in self.sample_names:
            self.process_one_pair(name)

        return self.out_fq_dir

    def set_sample_names(self):
        self.sample_names = pd.read_csv(self.sample_sheet, index_col=0).index.tolist()

    def set_out_fq_dir(self):
        self.out_fq_dir = f'{self.workdir}/trimmed_fastqs'
        os.makedirs(self.out_fq_dir, exist_ok=True)

    def process_one_pair(self, name: str):
        fq1 = f'{self.fq_dir}/{name}{self.fq1_suffix}'
        fq2 = f'{self.fq_dir}/{name}{self.fq2_suffix}'

        trimmed_fq1_gz, trimmed_fq2_gz = self.trim_galore(
            fq1=fq1,
            fq2=fq2,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime
        )

        for fq_gz, suffix in [
            (trimmed_fq1_gz, self.fq1_suffix),
            (trimmed_fq2_gz, self.fq2_suffix)
        ]:
            dst = f'{self.out_fq_dir}/{name}{suffix}'
            if suffix.endswith('.gz'):
                self.call(f'mv {fq_gz} {dst}')
            else:  # suffix is not .gz, so the dst file needs to be uncompressed
                self.call(f'gunzip {fq_gz} && mv {fq_gz[:-3]} {dst}')


class TrimGaloreSingleEnd(Processor):

    QUALITY = 20
    LENGTH = 20
    MAX_N = 0
    CUTADAPT_TOTAL_CORES = 2
    # According to the help message of trim_galore, 2 cores for cutadapt -> actually up to 9 cores

    fq: str
    clip_5_prime: int

    out_fq: str

    def main(self, fq: str, clip_5_prime: int) -> str:
        self.fq = fq
        self.clip_5_prime = clip_5_prime

        self.execute()
        self.move_fastqc_report()
        self.set_out_fq()

        return self.out_fq

    def execute(self):
        args = [
            'trim_galore',
            f'--quality {self.QUALITY}',
            '--phred33',
            f'--cores {self.CUTADAPT_TOTAL_CORES}',
            f'--fastqc_args "--threads {self.threads}"',
            '--illumina',
            f'--length {self.LENGTH}',
            f'--max_n {self.MAX_N}',
            '--trim-n',
            '--gzip',
            f'--output_dir {self.workdir}',
        ]

        if self.clip_5_prime > 0:
            args.append(f'--clip_R1 {self.clip_5_prime}')

        log = f'{self.outdir}/trim_galore.log'
        args += [
            self.fq,
            f'1>> {log} 2>> {log}'
        ]

        self.call(self.CMD_LINEBREAK.join(args))

    def move_fastqc_report(self):
        dstdir = f'{self.outdir}/fastqc'
        os.makedirs(dstdir, exist_ok=True)
        for suffix in [
            'fastqc.html',
            'fastqc.zip',
            'trimming_report.txt'
        ]:
            self.call(f'mv {self.workdir}/*{suffix} {dstdir}/')

    def set_out_fq(self):
        f = basename(self.fq)
        f = self.__strip_file_extension(f)
        self.out_fq = f'{self.workdir}/{f}_trimmed.fq.gz'

    def __strip_file_extension(self, f):
        for suffix in [
            '.fq',
            '.fq.gz',
            '.fastq',
            '.fastq.gz',
        ]:
            if f.endswith(suffix):
                f = f[:-len(suffix)]  # strip suffix
        return f


class BatchTrimGaloreSingleEnd(Processor):

    sample_sheet: str
    fq_dir: str
    fq_suffix: str
    clip_5_prime: int

    sample_names: List[str]
    out_fq_dir: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.trim_galore = TrimGaloreSingleEnd(self.settings).main

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq_suffix: str,
            clip_5_prime: int) -> str:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq_suffix = fq_suffix
        self.clip_5_prime = clip_5_prime

        self.set_sample_names()
        self.set_out_fq_dir()
        for name in self.sample_names:
            self.process_one_fq(name)

        return self.out_fq_dir

    def set_sample_names(self):
        self.sample_names = pd.read_csv(self.sample_sheet, index_col=0).index.tolist()

    def set_out_fq_dir(self):
        self.out_fq_dir = f'{self.workdir}/trimmed_fastqs'
        os.makedirs(self.out_fq_dir, exist_ok=True)

    def process_one_fq(self, name: str):
        fq = f'{self.fq_dir}/{name}{self.fq_suffix}'

        fq = self.trim_galore(fq=fq, clip_5_prime=self.clip_5_prime)

        os.rename(
            fq, f'{self.out_fq_dir}/{name}{self.fq_suffix}'
        )
