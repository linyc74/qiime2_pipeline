import pandas as pd
from .exporting import ExportTaxonomy
from .importing import ImportTaxonomy
from .template import Processor, Settings


class Taxonomy(Processor):

    CONFIDENCE = 0

    representative_seq_qza: str
    nb_classifier_qza: str

    forward_taxonomy_qza: str
    reverse_taxonomy_qza: str
    merged_taxonomy_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            representative_seq_qza: str,
            nb_classifier_qza: str) -> str:

        self.representative_seq_qza = representative_seq_qza
        self.nb_classifier_qza = nb_classifier_qza

        self.forward_classify()
        self.reverse_classify()
        self.merge_foward_reverse()

        return self.merged_taxonomy_qza

    def forward_classify(self):
        self.forward_taxonomy_qza = Classify(self.settings).main(
            representative_seq_qza=self.representative_seq_qza,
            nb_classifier_qza=self.nb_classifier_qza,
            read_orientation='same')

    def reverse_classify(self):
        self.reverse_taxonomy_qza = Classify(self.settings).main(
            representative_seq_qza=self.representative_seq_qza,
            nb_classifier_qza=self.nb_classifier_qza,
            read_orientation='reverse-complement')

    def merge_foward_reverse(self):
        self.merged_taxonomy_qza = MergeForwardReverseTaxonomy(self.settings).main(
            forward_taxonomy_qza=self.forward_taxonomy_qza,
            reverse_taxonomy_qza=self.reverse_taxonomy_qza)


class Classify(Processor):

    CONFIDENCE_CUTOFF = 0

    representative_seq_qza: str
    nb_classifier_qza: str
    read_orientation: str

    taxonomy_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            representative_seq_qza: str,
            nb_classifier_qza: str,
            read_orientation: str) -> str:

        self.representative_seq_qza = representative_seq_qza
        self.nb_classifier_qza = nb_classifier_qza
        self.read_orientation = read_orientation

        self.taxonomy_qza = f'{self.workdir}/taxonomy-{self.read_orientation}.qza'

        cmd = self.CMD_LINEBREAK.join([
            'qiime feature-classifier classify-sklearn',
            f'--i-classifier {self.nb_classifier_qza}',
            f'--i-reads {self.representative_seq_qza}',
            f'--p-read-orientation {self.read_orientation}',
            f'--p-confidence {self.CONFIDENCE_CUTOFF}',
            f'--p-n-jobs {self.threads}',
            f'--o-classification {self.taxonomy_qza}',
        ])
        self.call(cmd)

        return self.taxonomy_qza


class MergeForwardReverseTaxonomy(Processor):

    forward_taxonomy_qza: str
    reverse_taxonomy_qza: str

    f_df: pd.DataFrame
    r_df: pd.DataFrame
    df: pd.DataFrame

    merged_taxonomy_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(
            self,
            forward_taxonomy_qza: str,
            reverse_taxonomy_qza: str) -> str:

        self.forward_taxonomy_qza = forward_taxonomy_qza
        self.reverse_taxonomy_qza = reverse_taxonomy_qza

        self.read_taxonomy_qzas()
        self.rename_columns()
        self.merge_dfs()
        self.compare_by_confidence()
        self.save_as_qza()

        return self.merged_taxonomy_qza

    def read_taxonomy_qzas(self):
        tsv = ExportTaxonomy(self.settings).main(
            taxonomy_qza=self.forward_taxonomy_qza)
        self.f_df = pd.read_csv(tsv, sep='\t')
        tsv = ExportTaxonomy(self.settings).main(
            taxonomy_qza=self.reverse_taxonomy_qza)
        self.r_df = pd.read_csv(tsv, sep='\t')

    def rename_columns(self):
        self.f_df = self.f_df.rename(
            columns={
                'Taxon': 'Taxon (forward)',
                'Confidence': 'Confidence (forward)',
            }
        )
        self.r_df = self.r_df.rename(
            columns={
                'Taxon': 'Taxon (reverse)',
                'Confidence': 'Confidence (reverse)',
            }
        )

    def merge_dfs(self):
        self.df = self.f_df.merge(
            right=self.r_df,
            how='left',
            on='Feature ID'
        )
        assert len(self.df) == len(self.f_df)

    def compare_by_confidence(self):
        for i, row in self.df.iterrows():
            f_taxon = row['Taxon (forward)']
            r_taxon = row['Taxon (reverse)']
            f_confidence = row['Confidence (forward)']
            r_confidence = row['Confidence (reverse)']
            if f_confidence >= r_confidence:
                self.df.loc[i, 'Taxon'] = f_taxon
                self.df.loc[i, 'Confidence'] = f_confidence
            else:
                self.df.loc[i, 'Taxon'] = r_taxon
                self.df.loc[i, 'Confidence'] = r_confidence

        self.df.drop(
            columns=[
                'Taxon (forward)',
                'Taxon (reverse)',
                'Confidence (forward)',
                'Confidence (reverse)',
            ],
            inplace=True
        )

    def save_as_qza(self):
        tsv = f'{self.workdir}/taxonomy-merged.tsv'
        self.df.to_csv(tsv, sep='\t', index=False)
        self.merged_taxonomy_qza = ImportTaxonomy(self.settings).main(
            taxonomy_tsv=tsv)
