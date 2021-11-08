import os
import pandas as pd
from typing import List
from .tools import get_files
from .template import Processor, Settings


class ImportPairedEndFastq(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str

    manifest_tsv: str
    output_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str) -> str:

        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix

        self.manifest_tsv = f'{self.workdir}/manifest.tsv'
        self.output_qza = f'{self.workdir}/paired-end-demultiplexed.qza'

        self.write_manifest_tsv()
        self.import_with_manifest_tsv()

        return self.output_qza

    def write_manifest_tsv(self):
        WriteManifestTsv(self.settings).main(
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            manifest_tsv=self.manifest_tsv)

    def import_with_manifest_tsv(self):
        ImportWithManifestTsv(self.settings).main(
            manifest_tsv=self.manifest_tsv,
            output_qza=self.output_qza)


class WriteManifestTsv(Processor):

    fq_dir: str
    fq1_suffix: str
    fq2_suffix: str
    manifest_tsv: str

    sample_names: List[str]
    fq1_paths: List[str]
    fq2_paths: List[str]

    SAMPLE_ID_COLUMN = 'sample-id'
    FQ1_COLUMN = 'forward-absolute-filepath'
    FQ2_COLUMN = 'reverse-absolute-filepath'

    def __init__(self, settings: Settings):
        super().__init__(settings)

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

        self.set_sample_names()
        self.set_fq1_fq2_paths()
        self.assert_paths_exist()
        self.write_tsv()

    def set_sample_names(self):
        files = get_files(
            source=self.fq_dir,
            endswith=self.fq1_suffix,
            isfullpath=False)
        n_char = len(self.fq1_suffix)
        self.sample_names = [f[:-n_char] for f in files]

    def set_fq1_fq2_paths(self):
        self.fq1_paths = [
            f'{self.fq_dir}/{name}{self.fq1_suffix}'
            for name in self.sample_names
        ]
        self.fq2_paths = [
            f'{self.fq_dir}/{name}{self.fq2_suffix}'
            for name in self.sample_names
        ]

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


class ImportWithManifestTsv(Processor):

    TYPE = '\'SampleData[PairedEndSequencesWithQuality]\''
    INPUT_FORMAT = 'PairedEndFastqManifestPhred33V2'

    manifest_tsv: str
    output_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, manifest_tsv: str, output_qza: str):
        self.manifest_tsv = manifest_tsv
        self.output_qza = output_qza

        lines = [
            'qiime tools import',
            f'--type {self.TYPE}',
            f'--input-format {self.INPUT_FORMAT}',
            f'--input-path {self.manifest_tsv}',
            f'--output-path {self.output_qza}',
        ]
        self.call(' \\\n  '.join(lines))
