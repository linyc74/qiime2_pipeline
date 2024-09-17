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
    sample_sheet: str
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
            sample_sheet: str,
            skip_otu: bool) -> Tuple[str, str, str, str]:

        self.taxonomy_qza = taxonomy_qza
        self.feature_table_qza = feature_table_qza
        self.feature_sequence_qza = feature_sequence_qza
        self.sample_sheet = sample_sheet
        self.skip_otu = skip_otu

        self.decompress()
        self.set_feature_id_to_label()
        self.label_feature_sequence()
        self.label_feature_table()
        self.write_taxonomy_condifence_table()
        self.compress()

        return self.labeled_feature_table_tsv, self.labeled_feature_table_qza, \
            self.labeled_feature_sequence_fa, self.labeled_feature_sequence_qza

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
            feature_id_to_label=self.feature_id_to_label,
            sample_sheet=self.sample_sheet)

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
            taxon = ensure_one_space_after_semicolon(taxon)
            taxon = add_genus_prefix_to_unidentified_species(taxon)
            label = f'{self.label_prefix}{i + 1:04d}; {taxon}'
            self.output_dict[id_] = label

        return self.output_dict


def ensure_one_space_after_semicolon(taxon: str) -> str:
    items = taxon.split(';')
    items = [i.lstrip() for i in items]
    return '; '.join(items)


def add_genus_prefix_to_unidentified_species(taxon: str) -> str:
    if '; s__' not in taxon:
        return taxon  # species name is not available

    species = taxon.split('; s__')[-1]
    genus = taxon.split('; g__')[-1].split('; s__')[0]

    if genus.lower() in species.lower():
        return taxon

    prefix = taxon.split('; s__')[0]
    return f'{prefix}; s__{genus}_{species}'


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
    sample_sheet: str

    df: pd.DataFrame
    output_tsv: str

    def main(
            self,
            feature_table_tsv: str,
            feature_id_to_label: Dict[str, str],
            sample_sheet: str) -> str:

        self.feature_table_tsv = feature_table_tsv
        self.feature_id_to_label = feature_id_to_label
        self.sample_sheet = sample_sheet

        self.read_feature_table_tsv()
        self.reorder_sample_columns()
        self.convert_feature_id_to_label()
        self.save_output_tsv()

        return self.output_tsv

    def read_feature_table_tsv(self):
        self.df = pd.read_csv(
            self.feature_table_tsv,
            sep='\t',
            index_col=0)
        self.df.index.name = ''  # remove residual index name (e.g. "#OTU ID") created by biom

    def reorder_sample_columns(self):
        samples = pd.read_csv(self.sample_sheet, index_col=0).index
        self.df = self.df[samples]

    def convert_feature_id_to_label(self):
        self.df.index = [
            self.feature_id_to_label[idx] for idx in self.df.index
        ]

    def save_output_tsv(self):
        self.output_tsv = f'{self.outdir}/labeled-feature-table.tsv'
        self.df.to_csv(
            self.output_tsv,
            sep='\t',
            index=True)


class WriteTaxonomyCondifenceTable(Processor):

    taxonomy_tsv: str
    feature_id_to_label: Dict[str, str]

    df: pd.DataFrame
    confidence_or_consensus: str

    def main(
            self,
            taxonomy_tsv: str,
            feature_id_to_label: Dict[str, str]):

        self.taxonomy_tsv = taxonomy_tsv
        self.feature_id_to_label = feature_id_to_label

        self.read_taxonomy_tsv()
        self.set_confidence_or_consensus()
        self.label_and_process()
        self.save_output_tsv()

    def read_taxonomy_tsv(self):
        self.df = pd.read_csv(self.taxonomy_tsv, sep='\t')

    def set_confidence_or_consensus(self):
        if 'Confidence' in self.df.columns:
            self.confidence_or_consensus = 'Confidence'  # NB classifier
        elif 'Consensus' in self.df.columns:
            self.confidence_or_consensus = 'Consensus'  # VSEARCH classifier
        else:
            raise ValueError('No column named "Confidence" or "Consensus"')

    def label_and_process(self):
        self.df['Feature ID'] = [
            self.feature_id_to_label[id_] for id_ in self.df['Feature ID']
        ]
        self.df = self.df.rename(
            columns={'Feature ID': 'Feature Label'}
        )
        column = self.confidence_or_consensus
        self.df = self.df[['Feature Label', column]]

    def save_output_tsv(self):
        suffix = self.confidence_or_consensus.lower()
        output_tsv = f'{self.outdir}/taxonomy-{suffix}.tsv'
        self.df.to_csv(output_tsv, sep='\t', index=False)
