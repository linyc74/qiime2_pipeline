from typing import Tuple, Optional
from .template import Processor
from .pool import BatchPool
from .denoise import Dada2SingleEnd, Dada2PairedEnd
from .importing import ImportSingleEndFastq, ImportPairedEndFastq
from .trimming import BatchTrimGalorePairedEnd, BatchTrimGaloreSingleEnd


class GenerateASV(Processor):

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: Optional[str]
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
            paired_end_mode: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int,
            max_expected_error_bases: float) -> Tuple[str, str]:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.paired_end_mode = paired_end_mode
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime
        self.max_expected_error_bases = max_expected_error_bases

        if self.fq2_suffix is None:
            self.generate_asv_single_end()
        else:
            self.generate_asv_paired_end()

        return self.feature_table_qza, self.feature_sequence_qza

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
