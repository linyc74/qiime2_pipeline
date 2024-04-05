import os
import pandas as pd
from typing import List
from .tools import edit_fpath
from .template import Processor


class ImportFastq(Processor):

    sample_sheet: str
    fq_dir: str
    manifest_tsv: str
    output_qza: str


class ImportSingleEndFastq(ImportFastq):

    fq_suffix: str

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq_suffix: str) -> str:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq_suffix = fq_suffix

        self.write_manifest_tsv()
        self.import_with_manifest_tsv()

        return self.output_qza

    def write_manifest_tsv(self):
        self.manifest_tsv = f'{self.workdir}/manifest.tsv'
        WriteSingleEndManifestTsv(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq_suffix=self.fq_suffix,
            manifest_tsv=self.manifest_tsv)

    def import_with_manifest_tsv(self):
        self.output_qza = f'{self.workdir}/single-end-demultiplexed.qza'
        log = f'{self.outdir}/qiime-tools-import.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools import',
            f'--type \'SampleData[SequencesWithQuality]\'',
            f'--input-format SingleEndFastqManifestPhred33V2',
            f'--input-path {self.manifest_tsv}',
            f'--output-path {self.output_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)


class ImportPairedEndFastq(ImportFastq):

    sample_sheet: str
    fq1_suffix: str
    fq2_suffix: str

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str) -> str:

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix

        self.write_manifest_tsv()
        self.import_with_manifest_tsv()

        return self.output_qza

    def write_manifest_tsv(self):
        self.manifest_tsv = f'{self.workdir}/manifest.tsv'
        WritePairedEndManifestTsv(self.settings).main(
            sample_sheet=self.sample_sheet,
            fq_dir=self.fq_dir,
            fq1_suffix=self.fq1_suffix,
            fq2_suffix=self.fq2_suffix,
            manifest_tsv=self.manifest_tsv)

    def import_with_manifest_tsv(self):
        self.output_qza = f'{self.workdir}/paired-end-demultiplexed.qza'
        log = f'{self.outdir}/qiime-tools-import.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools import',
            f'--type \'SampleData[PairedEndSequencesWithQuality]\'',
            f'--input-format PairedEndFastqManifestPhred33V2',
            f'--input-path {self.manifest_tsv}',
            f'--output-path {self.output_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)


class WriteManifestTsv(Processor):

    # https://docs.qiime2.org/2021.11/tutorials/importing/
    SAMPLE_ID_COLUMN = 'sample-id'

    sample_sheet: str
    fq_dir: str
    manifest_tsv: str
    sample_names: List[str]

    def set_sample_names(self):
        self.sample_names = pd.read_csv(self.sample_sheet, index_col=0).index.tolist()


class WriteSingleEndManifestTsv(WriteManifestTsv):

    FQ_COLUMN = 'absolute-filepath'

    fq_suffix: str
    fq_paths: List[str]

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq_suffix: str,
            manifest_tsv: str):

        self.sample_sheet = sample_sheet
        self.fq_dir = os.path.abspath(fq_dir)  # qiime2 only accepts absolute path in the manifest tsv
        self.fq_suffix = fq_suffix
        self.manifest_tsv = manifest_tsv

        self.set_sample_names()
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
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: str,
            manifest_tsv: str):

        self.sample_sheet = sample_sheet
        self.fq_dir = os.path.abspath(fq_dir)  # qiime2 only accepts absolute path in the manifest tsv
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.manifest_tsv = manifest_tsv

        self.set_sample_names()
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


class ImportFeatureTable(Processor):

    feature_table_tsv: str
    biom: str
    qza: str

    def main(self, feature_table_tsv: str) -> str:
        self.feature_table_tsv = feature_table_tsv
        self.tsv_to_biom()
        self.biom_to_qza()
        return self.qza

    def tsv_to_biom(self):
        self.biom = edit_fpath(
            fpath=self.feature_table_tsv,
            old_suffix='.tsv',
            new_suffix='.biom',
            dstdir=self.workdir)
        log = f'{self.outdir}/biom-convert.log'
        cmd = self.CMD_LINEBREAK.join([
            'biom convert',
            '--to-hdf5',
            '--table-type="OTU table"',
            f'-i "{self.feature_table_tsv}"',
            f'-o "{self.biom}"',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)

    def biom_to_qza(self):
        self.qza = edit_fpath(
            fpath=self.biom,
            old_suffix='.biom',
            new_suffix='.qza',
            dstdir=self.workdir)
        log = f'{self.outdir}/qiime-tools-import.log'
        # https://biom-format.org/documentation/format_versions/biom-2.1.html
        # HDF5 is Biom file format version 2.1
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools import',
            f'--type \'FeatureTable[Frequency]\'',
            f'--input-format BIOMV210Format',
            f'--input-path "{self.biom}"',
            f'--output-path "{self.qza}"',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)


class ImportFeatureSequence(Processor):

    feature_sequence_fa: str
    qza: str

    def main(self, feature_sequence_fa: str) -> str:
        self.feature_sequence_fa = feature_sequence_fa
        self.fa_to_qza()
        return self.qza

    def fa_to_qza(self):
        self.qza = edit_fpath(
            fpath=self.feature_sequence_fa,
            old_suffix='.fa',
            new_suffix='.qza',
            dstdir=self.workdir)
        log = f'{self.outdir}/qiime-tools-import.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools import',
            f'--type FeatureData[Sequence]',
            f'--input-path "{self.feature_sequence_fa}"',
            f'--output-path "{self.qza}"',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)


class ImportTaxonomy(Processor):

    taxonomy_tsv: str

    qza: str

    def main(self, taxonomy_tsv: str) -> str:
        self.taxonomy_tsv = taxonomy_tsv
        self.tsv_to_qza()
        return self.qza

    def tsv_to_qza(self):
        self.qza = edit_fpath(
            fpath=self.taxonomy_tsv,
            old_suffix='.tsv',
            new_suffix='.qza',
            dstdir=self.workdir)
        log = f'{self.outdir}/qiime-tools-import.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools import',
            f'--type FeatureData[Taxonomy]',
            f'--input-path "{self.taxonomy_tsv}"',
            f'--output-path "{self.qza}"',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)
