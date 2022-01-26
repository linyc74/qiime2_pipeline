import pandas as pd
from typing import Tuple, Dict
from .template import Processor, Settings
from .fasta import FastaParser, FastaWriter
from .importing import ImportFeatureTable, ImportFeatureSequence
from .exporting import ExportFeatureTable, ExportFeatureSequence, ExportTaxonomy


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
