import os
import pandas as pd
from typing import Tuple, Optional, List
from .template import Processor
from .importing import ImportSingleEndFastq, ImportPairedEndFastq
from .trimming import BatchTrimGalorePairedEnd, BatchTrimGaloreSingleEnd


class GenerateASV(Processor):

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: Optional[str]
    pacbio: bool
    paired_end_mode: str
    clip_r1_5_prime: int
    clip_r2_5_prime: int
    max_expected_error_bases: float

    feature_table_qza: str
    feature_sequence_qza: str

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: Optional[str],
            pacbio: bool,
            paired_end_mode: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int,
            max_expected_error_bases: float) -> Tuple[str, str]:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.pacbio = pacbio
        self.paired_end_mode = paired_end_mode
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime
        self.max_expected_error_bases = max_expected_error_bases

        if self.pacbio:
            self.generate_asv_pacbio()

        else:  # illumina
            if self.fq2_suffix is None:
                self.generate_asv_single_end()
            else:
                self.generate_asv_paired_end()

        return self.feature_table_qza, self.feature_sequence_qza

    def generate_asv_pacbio(self):
        self.feature_table_qza, self.feature_sequence_qza = GenerateASVPacBio(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq_suffix=self.fq1_suffix,
            max_expected_error_bases=self.max_expected_error_bases)

    def generate_asv_single_end(self):
        self.feature_table_qza, self.feature_sequence_qza = GenerateASVSingleEnd(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq_suffix=self.fq1_suffix,
            clip_5_prime=self.clip_r1_5_prime,
            max_expected_error_bases=self.max_expected_error_bases)

    def generate_asv_paired_end(self):
        self.feature_table_qza, self.feature_sequence_qza = GenerateASVPairedEnd(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime,
            paired_end_mode=self.paired_end_mode,
            max_expected_error_bases=self.max_expected_error_bases)


class GenerateASVPacBio(Processor):

    sample_sheet: str
    fq_dir: str
    fq_suffix: str
    max_expected_error_bases: float

    feature_sequence_qza: str
    feature_table_qza: str

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq_suffix: str,
            max_expected_error_bases: float) -> Tuple[str, str]:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq_suffix = fq_suffix
        self.max_expected_error_bases = max_expected_error_bases

        # 1 importing
        single_end_seq_qza = ImportSingleEndFastq(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq_suffix=self.fq_suffix)

        # 2 denoise
        self.feature_table_qza, self.feature_sequence_qza = Dada2PacBio(self.settings).main(
            demultiplexed_seq_qza=single_end_seq_qza,
            max_expected_error_bases=self.max_expected_error_bases)

        return self.feature_table_qza, self.feature_sequence_qza


class GenerateASVSingleEnd(Processor):

    sample_sheet: str
    fq_dir: str
    fq_suffix: str
    clip_5_prime: int
    max_expected_error_bases: float

    feature_sequence_qza: str
    feature_table_qza: str

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq_suffix: str,
            clip_5_prime: int,
            max_expected_error_bases: float) -> Tuple[str, str]:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq_suffix = fq_suffix
        self.clip_5_prime = clip_5_prime
        self.max_expected_error_bases = max_expected_error_bases

        # 1 trimming
        trimmed_fq_dir, trimmed_fq_suffix = BatchTrimGaloreSingleEnd(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq_suffix=self.fq_suffix,
            clip_5_prime=self.clip_5_prime)

        # 2 importing
        single_end_seq_qza = ImportSingleEndFastq(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=trimmed_fq_dir,
            fq_suffix=trimmed_fq_suffix)

        # 3 denoise
        self.feature_table_qza, self.feature_sequence_qza = Dada2SingleEnd(self.settings).main(
            demultiplexed_seq_qza=single_end_seq_qza,
            max_expected_error_bases=self.max_expected_error_bases)

        return self.feature_table_qza, self.feature_sequence_qza


class GenerateASVPairedEnd(Processor):

    MERGE = 'merge'
    POOL = 'pool'

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str
    clip_r1_5_prime: int
    clip_r2_5_prime: int
    paired_end_mode: str
    max_expected_error_bases: float

    feature_sequence_qza: str
    feature_table_qza: str

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int,
            paired_end_mode: str,
            max_expected_error_bases: float) -> Tuple[str, str]:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime
        self.paired_end_mode = paired_end_mode
        self.max_expected_error_bases = max_expected_error_bases

        assert self.paired_end_mode in [self.MERGE, self.POOL], f'"{self.paired_end_mode}" is not a valid mode for GenerateASV'

        if self.paired_end_mode == self.MERGE:
            self.run_merge_mode()
        else:
            self.run_pool_mode()

        return self.feature_table_qza, self.feature_sequence_qza

    def run_merge_mode(self):
        # 1 trimming
        trimmed_fq_dir, trimmed_fq1_suffix, trimmed_fq2_suffix = BatchTrimGalorePairedEnd(
            self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime)

        # 2 importing
        paired_end_seq_qza = ImportPairedEndFastq(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=trimmed_fq_dir,
            fq1_suffix=trimmed_fq1_suffix,
            fq2_suffix=trimmed_fq2_suffix)

        # 3 denoise
        self.feature_table_qza, self.feature_sequence_qza = Dada2PairedEnd(self.settings).main(
            demultiplexed_seq_qza=paired_end_seq_qza,
            max_expected_error_bases=self.max_expected_error_bases)

    def run_pool_mode(self):
        # 1 trimming
        trimmed_fq_dir, trimmed_fq1_suffix, trimmed_fq2_suffix = BatchTrimGalorePairedEnd(
            self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime)

        # 2 pooling
        pooled_fq_dir, pooled_fq_suffix = BatchPool(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=trimmed_fq_dir,
            fq1_suffix=trimmed_fq1_suffix,
            fq2_suffix=trimmed_fq2_suffix)

        # 3 importing
        single_end_seq_qza = ImportSingleEndFastq(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=pooled_fq_dir,
            fq_suffix=pooled_fq_suffix)

        # 4 denoise
        self.feature_table_qza, self.feature_sequence_qza = Dada2SingleEnd(self.settings).main(
            demultiplexed_seq_qza=single_end_seq_qza,
            max_expected_error_bases=self.max_expected_error_bases)


#


class Dada2Base(Processor):

    TRIM_LEFT = 0  # number of 5' bases to be clipped
    TRUNCATE_LENGTH = 0  # min read length

    demultiplexed_seq_qza: str
    max_expected_error_bases: float

    feature_sequence_qza: str
    feature_table_qza: str
    denoising_stats_qza: str

    cmd: str

    def main(
            self,
            demultiplexed_seq_qza: str,
            max_expected_error_bases: float) -> Tuple[str, str]:

        self.demultiplexed_seq_qza = demultiplexed_seq_qza
        self.max_expected_error_bases = max_expected_error_bases

        self.set_output_paths()
        self.set_cmd()
        self.call(self.cmd)
        self.export_stats()

        return self.feature_table_qza, self.feature_sequence_qza

    def set_output_paths(self):
        self.feature_sequence_qza = f'{self.workdir}/dada2-feature-sequence.qza'
        self.feature_table_qza = f'{self.workdir}/dada2-feature-table.qza'
        self.denoising_stats_qza = f'{self.workdir}/dada2-stats.qza'

    def set_cmd(self):
        pass

    def export_stats(self):
        out = f'{self.workdir}/dada2'
        log = f'{self.outdir}/qiime-tools-export.log'

        cmd = self.CMD_LINEBREAK.join([
            'qiime tools export',
            f'--input-path {self.denoising_stats_qza}',
            f'--output-path {out}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)

        self.call(f'mv "{out}/stats.tsv" "{self.outdir}/dada2-stats.tsv"')


class Dada2SingleEnd(Dada2Base):

    def set_cmd(self):
        log = f'{self.outdir}/qiime-dada2-denoise-single.log'
        self.cmd = self.CMD_LINEBREAK.join([
            'qiime dada2 denoise-single',
            f'--i-demultiplexed-seqs {self.demultiplexed_seq_qza}',
            f'--p-trim-left {self.TRIM_LEFT}',
            f'--p-trunc-len {self.TRUNCATE_LENGTH}',
            f'--p-max-ee {self.max_expected_error_bases}',
            f'--p-n-threads {self.threads}',
            f'--o-representative-sequences {self.feature_sequence_qza}',
            f'--o-table {self.feature_table_qza}',
            f'--o-denoising-stats {self.denoising_stats_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])


class Dada2PairedEnd(Dada2Base):

    MIN_OVERLAP = 12

    def set_cmd(self):
        log = f'{self.outdir}/qiime-dada2-denoise-paired.log'
        self.cmd = self.CMD_LINEBREAK.join([
            'qiime dada2 denoise-paired',
            f'--i-demultiplexed-seqs {self.demultiplexed_seq_qza}',
            f'--p-trim-left-f {self.TRIM_LEFT}',
            f'--p-trim-left-r {self.TRIM_LEFT}',
            f'--p-trunc-len-f {self.TRUNCATE_LENGTH}',
            f'--p-trunc-len-r {self.TRUNCATE_LENGTH}',
            f'--p-max-ee-f {self.max_expected_error_bases}',
            f'--p-max-ee-r {self.max_expected_error_bases}',
            f'--p-min-overlap {self.MIN_OVERLAP}',
            f'--p-n-threads {self.threads}',
            f'--o-representative-sequences {self.feature_sequence_qza}',
            f'--o-table {self.feature_table_qza}',
            f'--o-denoising-stats {self.denoising_stats_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])


class Dada2PacBio(Dada2Base):
    """
    Reference:
      Callahan et al., High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide resolution (doi: 10.1093/nar/gkz569)
    """

    FORWARD_PRIMER = 'AGRGTTYGATYMTGGCTCAG'  # 27F
    REVERSE_PRIMER = 'RGYTACCTTGTTACGACTT'  # 1492R
    MIN_LENGTH = 1000
    MAX_LENGTH = 1600
    MIN_FOLD_PARENT_OVER_ABUNDANCE = 3.5

    def set_cmd(self):
        log = f'{self.outdir}/qiime-dada2-denoise-ccs.log'
        self.cmd = self.CMD_LINEBREAK.join([
            'qiime dada2 denoise-ccs',
            f'--i-demultiplexed-seqs {self.demultiplexed_seq_qza}',
            f'--p-front {self.FORWARD_PRIMER}',
            f'--p-adapter {self.REVERSE_PRIMER}',
            f'--p-max-ee {self.max_expected_error_bases}',
            f'--p-min-len {self.MIN_LENGTH}',
            f'--p-max-len {self.MAX_LENGTH}',
            f'--p-min-fold-parent-over-abundance {self.MIN_FOLD_PARENT_OVER_ABUNDANCE}',
            f'--p-n-threads {self.threads}',
            f'--o-representative-sequences {self.feature_sequence_qza}',
            f'--o-table {self.feature_table_qza}',
            f'--o-denoising-stats {self.denoising_stats_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])


#


class BatchPool(Processor):

    OUT_FQ_DIRNAME = 'pool_fastqs'
    OUT_FQ_SUFFIX = '.fastq.gz'

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
        self.make_out_fq_dir()
        for name in self.sample_names:
            self.process_one_pair(name=name)

        return self.out_fq_dir, self.OUT_FQ_SUFFIX

    def set_sample_names(self):
        self.sample_names = pd.read_csv(self.sample_sheet, index_col=0).index.tolist()

    def make_out_fq_dir(self):
        self.out_fq_dir = f'{self.workdir}/{self.OUT_FQ_DIRNAME}'
        os.makedirs(self.out_fq_dir, exist_ok=True)

    def process_one_pair(self, name: str):
        fq1 = f'{self.fq_dir}/{name}{self.fq1_suffix}'
        fq2 = f'{self.fq_dir}/{name}{self.fq2_suffix}'
        self.call(f'cat {fq1} {fq2} > {self.out_fq_dir}/{name}{self.OUT_FQ_SUFFIX}')
