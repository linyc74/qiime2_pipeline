from typing import Tuple
from .concat import BatchConcat
from .denoise import Dada2SingleEnd
from .trimming import BatchTrimGalore
from .template import Processor, Settings
from .importing import ImportSingleEndFastq


class GenerateASV(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str

    trimmed_fq_dir: str
    concat_fq_dir: str
    fq_suffix: str
    trimmed_reads_qza: str
    feature_sequence_qza: str
    feature_table_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str) -> Tuple[str, str]:

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix

        self.trimming()
        self.concat()
        self.importing()
        self.denoise()

        return self.feature_sequence_qza, self.feature_table_qza

    def trimming(self):
        self.trimmed_fq_dir = BatchTrimGalore(self.settings).main(
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def concat(self):
        self.concat_fq_dir, self.fq_suffix = BatchConcat(self.settings).main(
            fq_dir=self.trimmed_fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix)

    def importing(self):
        self.trimmed_reads_qza = ImportSingleEndFastq(self.settings).main(
            fq_dir=self.concat_fq_dir,
            fq_suffix=self.fq_suffix)

    def denoise(self):
        self.feature_sequence_qza, self.feature_table_qza = Dada2SingleEnd(self.settings).main(
            demultiplexed_seq_qza=self.trimmed_reads_qza)
