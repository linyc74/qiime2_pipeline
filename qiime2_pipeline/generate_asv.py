from typing import Tuple, Callable, Optional
from .template import Processor
from .concat import BatchConcat, BatchPool
from .denoise import Dada2SingleEnd, Dada2PairedEnd
from .importing import ImportSingleEndFastq, ImportPairedEndFastq
from .trimming import BatchTrimGalorePairedEnd, BatchTrimGaloreSingleEnd


class GenerateASV(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: Optional[str]
    paired_end_mode: str
    clip_r1_5_prime: int
    clip_r2_5_prime: int

    feature_table_qza: str
    feature_sequence_qza: str

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: Optional[str],
            paired_end_mode: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int) -> Tuple[str, str]:

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.paired_end_mode = paired_end_mode
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        if self.fq2_suffix is None:
            self.generate_asv_single_end()
        else:
            self.generate_asv_paired_end()

        return self.feature_table_qza, self.feature_sequence_qza

    def generate_asv_single_end(self):
        self.feature_table_qza, self.feature_sequence_qza = GenerateASVSingleEnd(self.settings).main(
            fq_dir=self.fq_dir,
            fq_suffix=self.fq1_suffix,
            clip_5_prime=self.clip_r1_5_prime)

    def generate_asv_paired_end(self):
        generate_asv = FactoryGenerateASVPairedEnd(self.settings).main(
            paired_end_mode=self.paired_end_mode)
        self.feature_table_qza, self.feature_sequence_qza = generate_asv(
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime)


class GenerateASVPairedEnd(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str
    clip_r1_5_prime: int
    clip_r2_5_prime: int

    trimmed_fq_dir: str
    feature_sequence_qza: str
    feature_table_qza: str

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int) -> Tuple[str, str]:

        return self.feature_table_qza, self.feature_sequence_qza

    def trimming(self):
        self.trimmed_fq_dir = BatchTrimGalorePairedEnd(self.settings).main(
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime)


class GenerateASVConcatPairedEnd(GenerateASVPairedEnd):

    concat_fq_dir: str
    fq_suffix: str
    single_end_seq_qza: str

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int) -> Tuple[str, str]:

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        self.trimming()
        self.concat()
        self.importing()
        self.denoise()

        return self.feature_table_qza, self.feature_sequence_qza

    def concat(self):
        self.concat_fq_dir, self.fq_suffix = BatchConcat(self.settings).main(
            fq_dir=self.trimmed_fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def importing(self):
        self.single_end_seq_qza = ImportSingleEndFastq(self.settings).main(
            fq_dir=self.concat_fq_dir,
            fq_suffix=self.fq_suffix)

    def denoise(self):
        self.feature_table_qza, self.feature_sequence_qza = Dada2SingleEnd(self.settings).main(
            demultiplexed_seq_qza=self.single_end_seq_qza)


class GenerateASVMergePairedEnd(GenerateASVPairedEnd):

    paired_end_seq_qza: str

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int) -> Tuple[str, str]:

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        self.trimming()
        self.importing()
        self.denoise()

        return self.feature_table_qza, self.feature_sequence_qza

    def importing(self):
        self.paired_end_seq_qza = ImportPairedEndFastq(self.settings).main(
            fq_dir=self.trimmed_fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def denoise(self):
        self.feature_table_qza, self.feature_sequence_qza = Dada2PairedEnd(self.settings).main(
            demultiplexed_seq_qza=self.paired_end_seq_qza)


class GenerateASVPoolPairedEnd(GenerateASVPairedEnd):

    pooled_fq_dir: str
    fq_suffix: str
    single_end_seq_qza: str

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int) -> Tuple[str, str]:

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        self.trimming()
        self.pool()
        self.importing()
        self.denoise()

        return self.feature_table_qza, self.feature_sequence_qza

    def pool(self):
        self.pooled_fq_dir, self.fq_suffix = BatchPool(self.settings).main(
            fq_dir=self.trimmed_fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def importing(self):
        self.single_end_seq_qza = ImportSingleEndFastq(self.settings).main(
            fq_dir=self.pooled_fq_dir,
            fq_suffix=self.fq_suffix)

    def denoise(self):
        self.feature_table_qza, self.feature_sequence_qza = Dada2SingleEnd(self.settings).main(
            demultiplexed_seq_qza=self.single_end_seq_qza)


class FactoryGenerateASVPairedEnd(Processor):

    MODE_TO_CLASS = {
        'concat': GenerateASVConcatPairedEnd,
        'merge': GenerateASVMergePairedEnd,
        'pool': GenerateASVPoolPairedEnd
    }

    def main(self, paired_end_mode: str) -> Callable:
        assert paired_end_mode in self.MODE_TO_CLASS.keys(), \
            f'"{paired_end_mode}" is not a valid mode for GenerateASV'

        _Class = self.MODE_TO_CLASS[paired_end_mode]
        return _Class(self.settings).main


class GenerateASVSingleEnd(Processor):

    fq_dir: str
    fq_suffix: str
    clip_5_prime: int

    trimmed_fq_dir: str
    single_end_seq_qza: str
    feature_sequence_qza: str
    feature_table_qza: str

    def main(
            self,
            fq_dir: str,
            fq_suffix: str,
            clip_5_prime: int) -> Tuple[str, str]:

        self.fq_dir = fq_dir
        self.fq_suffix = fq_suffix
        self.clip_5_prime = clip_5_prime

        self.trimming()
        self.importing()
        self.denoise()

        return self.feature_table_qza, self.feature_sequence_qza

    def trimming(self):
        self.trimmed_fq_dir = BatchTrimGaloreSingleEnd(self.settings).main(
            fq_dir=self.fq_dir,
            fq_suffix=self.fq_suffix,
            clip_5_prime=self.clip_5_prime)

    def importing(self):
        self.single_end_seq_qza = ImportSingleEndFastq(self.settings).main(
            fq_dir=self.trimmed_fq_dir,
            fq_suffix=self.fq_suffix)

    def denoise(self):
        self.feature_table_qza, self.feature_sequence_qza = Dada2SingleEnd(self.settings).main(
            demultiplexed_seq_qza=self.single_end_seq_qza)
