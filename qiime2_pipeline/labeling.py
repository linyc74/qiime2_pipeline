import pandas as pd
from typing import Tuple, Dict
from .template import Processor
from .fasta import FastaParser, FastaWriter
from .importing import ImportFeatureTable, ImportFeatureSequence
from .exporting import ExportFeatureTable, ExportFeatureSequence, ExportTaxonomy


class FeatureLabeling(Processor):

    ASV_PREFIX = 'ASV_'
    OTU_PREFIX = 'OTU_'

    taxonomy_qza: str
    feature_table_qza: str
    feature_sequence_qza: str
    skip_otu: bool

    taxonomy_tsv: str
    feature_table_tsv: str
    feature_sequence_fa: str

    feature_id_to_label: Dict[str, str]

    labeled_feature_table_tsv: str
    labeled_feature_sequence_fa: str

    labeled_feature_table_qza: str
    labeled_feature_sequence_qza: str

    def main(
            self,
            taxonomy_qza: str,
            feature_table_qza: str,
            feature_sequence_qza: str,
            skip_otu: bool) -> Tuple[str, str]:

        self.taxonomy_qza = taxonomy_qza
        self.feature_table_qza = feature_table_qza
        self.feature_sequence_qza = feature_sequence_qza
        self.skip_otu = skip_otu

        self.decompress()
        self.set_feature_id_to_label()
        self.label_feature_sequence()
        self.label_feature_table()
        self.write_taxonomy_condifence_table()
        self.compress()

        return self.labeled_feature_table_qza, self.labeled_feature_sequence_qza

    def decompress(self):
        self.taxonomy_tsv = ExportTaxonomy(self.settings).main(
            taxonomy_qza=self.taxonomy_qza)

        self.feature_table_tsv = ExportFeatureTable(self.settings).main(
            feature_table_qza=self.feature_table_qza)

        self.feature_sequence_fa = ExportFeatureSequence(self.settings).main(
            feature_sequence_qza=self.feature_sequence_qza)

    def set_feature_id_to_label(self):
        label_prefix = self.ASV_PREFIX if self.skip_otu else self.OTU_PREFIX
        self.feature_id_to_label = GetFeatureIDToLabelDict(self.settings).main(
            taxonomy_tsv=self.taxonomy_tsv,
            label_prefix=label_prefix)

    def label_feature_sequence(self):
        self.labeled_feature_sequence_fa = LabelFeatureSequence(self.settings).main(
            feature_sequence_fa=self.feature_sequence_fa,
            feature_id_to_label=self.feature_id_to_label)

    def label_feature_table(self):
        self.labeled_feature_table_tsv = LabelFeatureTable(self.settings).main(
            feature_table_tsv=self.feature_table_tsv,
            feature_id_to_label=self.feature_id_to_label)

    def write_taxonomy_condifence_table(self):
        WriteTaxonomyCondifenceTable(self.settings).main(
            taxonomy_tsv=self.taxonomy_tsv,
            feature_id_to_label=self.feature_id_to_label)

    def compress(self):
        self.labeled_feature_table_qza = ImportFeatureTable(self.settings).main(
            feature_table_tsv=self.labeled_feature_table_tsv)

        self.labeled_feature_sequence_qza = ImportFeatureSequence(self.settings).main(
            feature_sequence_fa=self.labeled_feature_sequence_fa)


class GetFeatureIDToLabelDict(Processor):

    taxonomy_tsv: str
    label_prefix: str

    df: pd.DataFrame
    output_dict: Dict[str, str]

    def main(
            self,
            taxonomy_tsv: str,
            label_prefix: str) -> Dict[str, str]:

        self.taxonomy_tsv = taxonomy_tsv
        self.label_prefix = label_prefix

        self.df = pd.read_csv(self.taxonomy_tsv, sep='\t')

        self.output_dict = {}
        for i, row in self.df.iterrows():
            id_ = row['Feature ID']
            taxon = row['Taxon']
            label = f'{self.label_prefix}{i + 1:04d}; {taxon}'
            self.output_dict[id_] = label

        return self.output_dict


class LabelFeatureSequence(Processor):

    feature_sequence_fa: str
    feature_id_to_label: Dict[str, str]

    output_fa: str

    def main(
            self,
            feature_sequence_fa: str,
            feature_id_to_label: Dict[str, str]) -> str:

        self.feature_sequence_fa = feature_sequence_fa
        self.feature_id_to_label = feature_id_to_label

        self.output_fa = f'{self.outdir}/labeled-feature-sequence.fa'

        with FastaParser(self.feature_sequence_fa) as parser:
            with FastaWriter(self.output_fa) as writer:
                for id_, seq in parser:
                    label = self.feature_id_to_label[id_]
                    writer.write(label, seq)

        return self.output_fa


class LabelFeatureTable(Processor):

    feature_table_tsv: str
    feature_id_to_label: Dict[str, str]

    df: pd.DataFrame
    output_tsv: str

    def main(
            self,
            feature_table_tsv: str,
            feature_id_to_label: Dict[str, str]) -> str:

        self.feature_table_tsv = feature_table_tsv
        self.feature_id_to_label = feature_id_to_label

        self.read_feature_table_tsv()
        self.convert_feature_id_to_label()
        self.save_output_tsv()

        return self.output_tsv

    def read_feature_table_tsv(self):
        self.df = pd.read_csv(
            self.feature_table_tsv,
            sep='\t',
            skiprows=1  # exclude 1st line from qza (# Constructed from biom file)
        )

    def convert_feature_id_to_label(self):
        self.df['#OTU ID'] = [
            self.feature_id_to_label[id_]
            for id_ in self.df['#OTU ID']
        ]
        self.df.rename(
            columns={'#OTU ID': ''},
            inplace=True
        )

    def save_output_tsv(self):
        self.output_tsv = f'{self.outdir}/labeled-feature-table.tsv'
        self.df.to_csv(
            self.output_tsv,
            sep='\t',
            index=False)


class WriteTaxonomyCondifenceTable(Processor):

    taxonomy_tsv: str
    feature_id_to_label: Dict[str, str]

    df: pd.DataFrame

    def main(
            self,
            taxonomy_tsv: str,
            feature_id_to_label: Dict[str, str]):

        self.taxonomy_tsv = taxonomy_tsv
        self.feature_id_to_label = feature_id_to_label

        self.read_taxonomy_tsv()
        self.label_and_process()
        self.save_output_tsv()

    def read_taxonomy_tsv(self):
        self.df = pd.read_csv(self.taxonomy_tsv, sep='\t')

    def label_and_process(self):
        self.df['Feature ID'] = [
            self.feature_id_to_label[id_]
            for id_ in self.df['Feature ID']
        ]
        self.df = self.df.rename(
            columns={'Feature ID': 'Feature Label'}
        )
        self.df = self.df[['Feature Label', 'Confidence']]

    def save_output_tsv(self):
        output_tsv = f'{self.outdir}/taxonomy-condifence.tsv'
        self.df.to_csv(output_tsv, sep='\t', index=False)
