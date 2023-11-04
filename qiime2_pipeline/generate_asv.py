from typing import Tuple, Optional
from .template import Processor
from .pool import BatchPool
from .denoise import Dada2SingleEnd, Dada2PairedEnd
from .importing import ImportSingleEndFastq, ImportPairedEndFastq
from .trimming import BatchTrimGalorePairedEnd, BatchTrimGaloreSingleEnd


class GenerateASV(Processor):

    MERGE = 'merge'
    POOL = 'pool'

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
        assert self.paired_end_mode in [self.MERGE, self.POOL], \
            f'"{self.paired_end_mode}" is not a valid mode for GenerateASV'

        if self.paired_end_mode == self.MERGE:
            method = GenerateASVMergePairedEnd(self.settings).main
        else:
            method = GenerateASVPoolPairedEnd(self.settings).main

        self.feature_table_qza, self.feature_sequence_qza = method(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime,
            max_expected_error_bases=self.max_expected_error_bases)


class GenerateASVSingleEnd(Processor):

    sample_sheet: str
    fq_dir: str
    fq_suffix: str
    clip_5_prime: int
    max_expected_error_bases: float

    trimmed_fq_dir: str
    single_end_seq_qza: str
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

        self.trimming()
        self.importing()
        self.denoise()

        return self.feature_table_qza, self.feature_sequence_qza

    def trimming(self):
        self.trimmed_fq_dir = BatchTrimGaloreSingleEnd(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq_suffix=self.fq_suffix,
            clip_5_prime=self.clip_5_prime)

    def importing(self):
        self.single_end_seq_qza = ImportSingleEndFastq(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.trimmed_fq_dir,
            fq_suffix=self.fq_suffix)

    def denoise(self):
        self.feature_table_qza, self.feature_sequence_qza = Dada2SingleEnd(self.settings).main(
            demultiplexed_seq_qza=self.single_end_seq_qza,
            max_expected_error_bases=self.max_expected_error_bases)


class GenerateASVPairedEnd(Processor):

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str
    clip_r1_5_prime: int
    clip_r2_5_prime: int
    max_expected_error_bases: float

    trimmed_fq_dir: str
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
            max_expected_error_bases: float) -> Tuple[str, str]:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime
        self.max_expected_error_bases = max_expected_error_bases

        self.workflow()

        return self.feature_table_qza, self.feature_sequence_qza

    def workflow(self):
        pass

    def trimming(self):
        self.trimmed_fq_dir = BatchTrimGalorePairedEnd(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime)


class GenerateASVMergePairedEnd(GenerateASVPairedEnd):

    paired_end_seq_qza: str

    def workflow(self):

        self.trimming()
        self.importing()
        self.denoise()

        return self.feature_table_qza, self.feature_sequence_qza

    def importing(self):
        self.paired_end_seq_qza = ImportPairedEndFastq(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.trimmed_fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def denoise(self):
        self.feature_table_qza, self.feature_sequence_qza = Dada2PairedEnd(self.settings).main(
            demultiplexed_seq_qza=self.paired_end_seq_qza,
            max_expected_error_bases=self.max_expected_error_bases)


class GenerateASVPoolPairedEnd(GenerateASVPairedEnd):

    pooled_fq_dir: str
    fq_suffix: str
    single_end_seq_qza: str

    def workflow(self):

        self.trimming()
        self.pool()
        self.importing()
        self.denoise()

        return self.feature_table_qza, self.feature_sequence_qza

    def pool(self):
        self.pooled_fq_dir, self.fq_suffix = BatchPool(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.trimmed_fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def importing(self):
        self.single_end_seq_qza = ImportSingleEndFastq(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.pooled_fq_dir,
            fq_suffix=self.fq_suffix)

    def denoise(self):
        self.feature_table_qza, self.feature_sequence_qza = Dada2SingleEnd(self.settings).main(
            demultiplexed_seq_qza=self.single_end_seq_qza,
            max_expected_error_bases=self.max_expected_error_bases)
