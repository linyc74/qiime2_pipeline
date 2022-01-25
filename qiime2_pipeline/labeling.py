import pandas as pd
from os.path import basename
from typing import Tuple, Dict
from .template import Processor, Settings
from .fasta import FastaParser, FastaWriter


class FeatureLabeling(Processor):

    taxonomy_qza: str
    feature_table_qza: str
    feature_sequence_qza: str

    taxonomy_tsv: str
    feature_table_tsv: str
    feature_sequence_fa: str

    feature_id_to_label: Dict[str, str]
    labeled_feature_table_tsv: str
    labeled_feature_sequence_fa: str

    labeled_feature_table_qza: str
    labeled_feature_sequence_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            taxonomy_qza: str,
            feature_table_qza: str,
            feature_sequence_qza: str) -> Tuple[str, str]:

        self.taxonomy_qza = taxonomy_qza
        self.feature_table_qza = feature_table_qza
        self.feature_sequence_qza = feature_sequence_qza

        self.decompress_input_qza()
        self.set_feature_id_to_label()
        self.label_feature_sequence()
        self.label_feature_table()
        self.write_taxon_condifence_table()
        self.compress_output_qza()

        return self.labeled_feature_table_qza, self.labeled_feature_sequence_qza

    def decompress_input_qza(self):
        self.taxonomy_tsv = ExportTaxonomy(self.settings).main(
            taxonomy_qza=self.taxonomy_qza)

        self.feature_table_tsv = ExportFeatureTable(self.settings).main(
            feature_table_qza=self.feature_table_qza)

        self.feature_sequence_fa = ExportFeatureSequence(self.settings).main(
            feature_sequence_qza=self.feature_sequence_qza)

    def set_feature_id_to_label(self):
        self.feature_id_to_label = GetFeatureIDToLabelDict(self.settings).main(
            taxonomy_tsv=self.taxonomy_tsv)

    def label_feature_sequence(self):
        self.labeled_feature_sequence_fa = f'{self.outdir}/labeled-feature-sequence.fa'

        with FastaParser(self.feature_sequence_fa) as parser:
            with FastaWriter(self.labeled_feature_sequence_fa) as writer:
                for id_, seq in parser:
                    label = self.feature_id_to_label[id_]
                    writer.write(label, seq)

    def label_feature_table(self):
        df = pd.read_csv(
            self.feature_table_tsv,
            sep='\t',
            skiprows=1
        )

        df['#OTU ID'] = [
            self.feature_id_to_label[id_]
            for id_ in df['#OTU ID']
        ]

        self.labeled_feature_table_tsv = f'{self.outdir}/labeled-feature-table.tsv'

        df.rename(
            columns={'#OTU ID': 'Feature Label'}
        ).to_csv(
            self.labeled_feature_table_tsv,
            sep='\t',
            index=False
        )

    def write_taxon_condifence_table(self):
        pass

    def compress_output_qza(self):
        self.labeled_feature_table_qza = ImportFeatureTable(self.settings).main(
            feature_table_tsv=self.labeled_feature_table_tsv)

        self.labeled_feature_sequence_qza = ImportFeatureSequence(self.settings).main(
            feature_sequence_fa=self.labeled_feature_sequence_fa)


class ExportFeatureTable(Processor):

    feature_table_qza: str
    output_tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_table_qza: str) -> str:
        self.feature_table_qza = feature_table_qza

        self.qza_to_biom()
        self.biom_to_tsv()

        return self.output_tsv

    def qza_to_biom(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.feature_table_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def biom_to_tsv(self):
        fname = basename(self.feature_table_qza)[:-len('.qza')] + '.tsv'
        self.output_tsv = f'{self.outdir}/{fname}'
        lines = [
            'biom convert --to-tsv',
            f'-i {self.workdir}/feature-table.biom',
            f'-o {self.output_tsv}'
        ]
        self.call(' \\\n  '.join(lines))


class ExportFeatureSequence(Processor):

    feature_sequence_qza: str
    output_fa: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_sequence_qza: str) -> str:

        self.feature_sequence_qza = feature_sequence_qza

        self.qza_to_fa()
        self.move_fa()

        return self.output_fa

    def qza_to_fa(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.feature_sequence_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def move_fa(self):
        fname = basename(self.feature_sequence_qza)[:-len('.qza')] + '.fa'
        self.output_fa = f'{self.outdir}/{fname}'
        self.call(f'mv {self.workdir}/dna-sequences.fasta {self.output_fa}')


class ExportTaxonomy(Processor):

    taxonomy_qza: str
    output_tsv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, taxonomy_qza: str) -> str:

        self.taxonomy_qza = taxonomy_qza

        self.set_output_tsv()
        self.qza_to_tsv()
        self.move_tsv()

        return self.output_tsv

    def qza_to_tsv(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.taxonomy_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def set_output_tsv(self):
        fname = basename(self.taxonomy_qza)[:-len('.qza')] + '.tsv'
        self.output_tsv = f'{self.workdir}/{fname}'

    def move_tsv(self):
        self.call(f'mv {self.workdir}/taxonomy.tsv {self.output_tsv}')


class GetFeatureIDToLabelDict(Processor):

    LABEL_PREFIX = 'ASV_'

    taxonomy_tsv: str

    df: pd.DataFrame
    output_dict: Dict[str, str]

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, taxonomy_tsv: str) -> Dict[str, str]:
        self.taxonomy_tsv = taxonomy_tsv

        self.df = pd.read_csv(self.taxonomy_tsv, sep='\t')

        self.output_dict = {}
        for i, row in self.df.iterrows():
            id_ = row['Feature ID']
            taxon = row['Taxon']
            label = f'{self.LABEL_PREFIX}{i+1:04d}; {taxon}'
            self.output_dict[id_] = label

        return self.output_dict


class ImportFeatureTable(Processor):

    feature_table_tsv: str
    biom: str
    output_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_table_tsv: str) -> str:
        self.feature_table_tsv = feature_table_tsv

        self.tsv_to_biom()
        self.biom_to_qza()

        return self.output_qza

    def tsv_to_biom(self):
        self.biom = edit_fpath(
            fpath=self.feature_table_tsv,
            old_suffix='.tsv',
            new_suffix='.biom',
            dstdir=self.workdir)

        cmd = self.CMD_LINEBREAK.join([
            'biom convert',
            '--to-hdf5',
            '--table-type="OTU table"',
            f'-i {self.feature_table_tsv}',
            f'-o {self.biom}'
        ])
        self.call(cmd)

    def biom_to_qza(self):
        self.output_qza = edit_fpath(
            fpath=self.biom,
            old_suffix='.biom',
            new_suffix='.qza',
            dstdir=self.workdir)

        # https://biom-format.org/documentation/format_versions/biom-2.1.html
        # HDF5 is the Biom file format version 2.1
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools import',
            f'--type \'FeatureTable[Frequency]\'',
            f'--input-format BIOMV210Format',
            f'--input-path {self.biom}',
            f'--output-path {self.output_qza}',
        ])
        self.call(cmd)


class ImportFeatureSequence(Processor):

    feature_sequence_fa: str
    feature_sequence_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, feature_sequence_fa: str) -> str:
        self.feature_sequence_fa = feature_sequence_fa

        self.set_feature_sequence_qza()
        self.execute()

        return self.feature_sequence_qza

    def execute(self):
        cmd = self.CMD_LINEBREAK.join([
            'qiime tools import',
            f'--type FeatureData[Sequence]',
            f'--input-path {self.feature_sequence_fa}',
            f'--output-path {self.feature_sequence_qza}',
        ])
        self.call(cmd)

    def set_feature_sequence_qza(self):
        self.feature_sequence_qza = edit_fpath(
            fpath=self.feature_sequence_fa,
            old_suffix='.fa',
            new_suffix='.qza',
            dstdir=self.workdir)


def edit_fpath(
        fpath: str,
        old_suffix: str,
        new_suffix: str,
        dstdir: str) -> str:

    fname = basename(fpath)[:-len(old_suffix)] + new_suffix
    return f'{dstdir}/{fname}'
