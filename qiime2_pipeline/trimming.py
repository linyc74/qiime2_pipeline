import os
from typing import Tuple
from os.path import basename
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


class TrimGalore(Processor):

    QUALITY = 20
    LENGTH = 20
    MAX_N = 0
    CUTADAPT_TOTAL_CORES = 2
    # According to the help message of trim_galore, 2 cores for cutadapt -> actually up to 9 cores

    fq1: str
    fq2: str

    out_fq1: str
    out_fq2: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, fq1: str, fq2: str) -> Tuple[str, str]:
        self.fq1 = fq1
        self.fq2 = fq2

        self.execute()
        self.move_fastqc_report()
        self.set_out_fq1()
        self.set_out_fq2()

        return self.out_fq1, self.out_fq2

    def execute(self):
        log = f'{self.outdir}/trim_galore.log'
        cmd = self.CMD_LINEBREAK.join([
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
            f'--output_dir {self.workdir}',
            self.fq1,
            self.fq2,
            f'1> {log} 2> {log}'
        ])
        self.call(cmd)

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
