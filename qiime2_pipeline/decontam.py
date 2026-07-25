import pandas as pd
from typing import Tuple
from .template import Processor


class Decontam(Processor):

    feature_table_qza: str
    feature_sequence_qza: str
    sample_sheet: str
    dna_concentration_column: str
    decontam_threshold: float

    concentration_tsv: str

    decontam_scores_qza: str
    decontam_table_qza: str
    decontam_sequence_qza: str

    def main(
            self,
            feature_table_qza: str,
            feature_sequence_qza: str,
            sample_sheet: str,
            dna_concentration_column: str,
            decontam_threshold: float) -> Tuple[str, str]:

        self.feature_table_qza = feature_table_qza
        self.feature_sequence_qza = feature_sequence_qza
        self.sample_sheet = sample_sheet
        self.dna_concentration_column = dna_concentration_column
        self.decontam_threshold = decontam_threshold

        self.write_concentration_tsv()
        self.identify_scores()
        self.visualize_score_histogram()
        self.remove_from_feature_table()
        self.remove_from_feature_sequences()

        return self.decontam_table_qza, self.decontam_sequence_qza

    def write_concentration_tsv(self):
        if self.sample_sheet.endswith('.csv'):
            df = pd.read_csv(self.sample_sheet, index_col=0)
        else:
            df = pd.read_csv(self.sample_sheet, sep='\t', index_col=0)

        df = df[[self.dna_concentration_column]]  # keep only one column needed

        # these two column names are required by qiime2 decontam-identify
        df.index.name = 'sample-id'
        df = df.rename(columns={self.dna_concentration_column: 'dna-concentration'})

        self.concentration_tsv = f'{self.workdir}/dna-concentration.tsv'
        df.to_csv(self.concentration_tsv, sep='\t', index=True)

    def identify_scores(self):
        self.decontam_scores_qza = f'{self.workdir}/decontam-scores.qza'
        log = f'{self.outdir}/qiime-quality-control-decontam-identify.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime quality-control decontam-identify',
            f'--i-table {self.feature_table_qza}',
            f'--m-metadata-file {self.concentration_tsv}',
            f'--p-method frequency',
            f'--p-freq-concentration-column "dna-concentration"',
            f'--o-decontam-scores {self.decontam_scores_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)

    def visualize_score_histogram(self):
        output_qzv = f'{self.outdir}/decontam-scores.qzv'
        log = f'{self.outdir}/qiime-quality-control-decontam-score-viz.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime quality-control decontam-score-viz',
            f'--i-decontam-scores {self.decontam_scores_qza}',
            f'--i-table {self.feature_table_qza}',
            f'--p-threshold {self.decontam_threshold}',
            f'--o-visualization {output_qzv}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)

    def remove_from_feature_table(self):
        self.decontam_table_qza = f'{self.workdir}/decontam-feature-table.qza'
        log = f'{self.outdir}/qiime-quality-control-decontam-remove.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime quality-control decontam-remove',
            f'--i-decontam-scores {self.decontam_scores_qza}',
            f'--i-table {self.feature_table_qza}',
            f'--p-threshold {self.decontam_threshold}',
            f'--o-filtered-table {self.decontam_table_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)
    
    def remove_from_feature_sequences(self):
        self.decontam_sequence_qza = f'{self.workdir}/decontam-feature-sequence.qza'
        log = f'{self.outdir}/qiime-feature-table-filter-seqs.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime feature-table filter-seqs',
            f'--i-data {self.feature_sequence_qza}',
            f'--i-table {self.feature_table_qza}',
            f'--o-filtered-data {self.decontam_sequence_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)
