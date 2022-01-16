import os
import pandas as pd
from typing import List
from .tools import get_files
from .template import Processor, Settings


class ImportFastq(Processor):

    fq_dir: str
    manifest_tsv: str
    output_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)


class ImportSingleEndFastq(ImportFastq):

    fq_suffix: str

    def main(
            self,
            fq_dir: str,
            fq_suffix: str) -> str:

        self.fq_dir = fq_dir
        self.fq_suffix = fq_suffix

        self.write_manifest_tsv()
        self.import_with_manifest_tsv()

        return self.output_qza

    def write_manifest_tsv(self):
        self.manifest_tsv = f'{self.workdir}/manifest.tsv'
        WriteSingleEndManifestTsv(self.settings).main(
            fq_dir=self.fq_dir,
            fq_suffix=self.fq_suffix,
            manifest_tsv=self.manifest_tsv)

    def import_with_manifest_tsv(self):
        self.output_qza = f'{self.workdir}/single-end-demultiplexed.qza'
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools import',
            f'--type \'SampleData[SequencesWithQuality]\'',
            f'--input-format SingleEndFastqManifestPhred33V2',
            f'--input-path {self.manifest_tsv}',
            f'--output-path {self.output_qza}',
        ])
        self.call(cmd)


class ImportPairedEndFastq(ImportFastq):

    fq1_suffix: str
    fq2_suffix: str

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str) -> str:

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix

        self.write_manifest_tsv()
        self.import_with_manifest_tsv()

        return self.output_qza

    def write_manifest_tsv(self):
        self.manifest_tsv = f'{self.workdir}/manifest.tsv'
        WritePairedEndManifestTsv(self.settings).main(
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            manifest_tsv=self.manifest_tsv)

    def import_with_manifest_tsv(self):
        self.output_qza = f'{self.workdir}/paired-end-demultiplexed.qza'
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools import',
            f'--type \'SampleData[PairedEndSequencesWithQuality]\'',
            f'--input-format PairedEndFastqManifestPhred33V2',
            f'--input-path {self.manifest_tsv}',
            f'--output-path {self.output_qza}',
        ])
        self.call(cmd)


class WriteManifestTsv(Processor):

    # https://docs.qiime2.org/2021.11/tutorials/importing/
    SAMPLE_ID_COLUMN = 'sample-id'

    fq_dir: str
    manifest_tsv: str
    sample_names: List[str]

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def set_sample_names(self, fq_suffix: str):
        files = get_files(
            source=self.fq_dir,
            endswith=fq_suffix,
            isfullpath=False)
        n_char = len(fq_suffix)
        self.sample_names = [f[:-n_char] for f in files]


class WriteSingleEndManifestTsv(WriteManifestTsv):

    FQ_COLUMN = 'absolute-filepath'

    fq_suffix: str
    fq_paths: List[str]

    def main(
            self,
            fq_dir: str,
            fq_suffix: str,
            manifest_tsv: str):

        self.fq_dir = os.path.abspath(fq_dir)  # qiime2 only accepts absolute path in the manifest tsv
        self.fq_suffix = fq_suffix
        self.manifest_tsv = manifest_tsv

        self.set_sample_names(fq_suffix=self.fq_suffix)
        self.set_fq_paths()
        self.write_tsv()

    def set_fq_paths(self):
        self.fq_paths = [
            f'{self.fq_dir}/{name}{self.fq_suffix}'
            for name in self.sample_names
        ]
        self.assert_paths_exist()

    def assert_paths_exist(self):
        for p in self.fq_paths:
            assert os.path.exists(p), f'{p} does not exist'

    def write_tsv(self):
        df = pd.DataFrame(data={
            self.SAMPLE_ID_COLUMN: self.sample_names,
            self.FQ_COLUMN: self.fq_paths,
        })
        df.to_csv(self.manifest_tsv, index=False, sep='\t')


class WritePairedEndManifestTsv(WriteManifestTsv):

    FQ1_COLUMN = 'forward-absolute-filepath'
    FQ2_COLUMN = 'reverse-absolute-filepath'

    fq1_suffix: str
    fq2_suffix: str

    fq1_paths: List[str]
    fq2_paths: List[str]

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            manifest_tsv: str):

        self.fq_dir = os.path.abspath(fq_dir)  # qiime2 only accepts absolute path in the manifest tsv
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.manifest_tsv = manifest_tsv

        self.set_sample_names(fq_suffix=self.fq1_suffix)
        self.set_fq1_fq2_paths()
        self.write_tsv()

    def set_fq1_fq2_paths(self):
        self.fq1_paths = [
            f'{self.fq_dir}/{name}{self.fq1_suffix}'
            for name in self.sample_names
        ]
        self.fq2_paths = [
            f'{self.fq_dir}/{name}{self.fq2_suffix}'
            for name in self.sample_names
        ]
        self.assert_paths_exist()

    def assert_paths_exist(self):
        for paths in [self.fq1_paths, self.fq2_paths]:
            for p in paths:
                assert os.path.exists(p), f'{p} does not exist'

    def write_tsv(self):
        df = pd.DataFrame(data={
            self.SAMPLE_ID_COLUMN: self.sample_names,
            self.FQ1_COLUMN: self.fq1_paths,
            self.FQ2_COLUMN: self.fq2_paths
        })
        df.to_csv(self.manifest_tsv, index=False, sep='\t')
