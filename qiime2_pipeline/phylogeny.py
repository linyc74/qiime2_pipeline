from typing import Tuple
from .template import Processor, Settings


class MafftFasttree(Processor):

    seq_qza: str

    aligned_seq_qza: str
    masked_aligned_seq_qza: str
    unrooted_tree_qza: str
    rooted_tree_qza: str
    output_paths: Tuple[str, str, str, str]

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, seq_qza: str) -> Tuple[str, str, str, str]:
        self.seq_qza = seq_qza

        self.set_output_paths()
        self.execute()

        return self.output_paths

    def set_output_paths(self):
        self.aligned_seq_qza = f'{self.workdir}/mafft-aligned-sequences.qza'
        self.masked_aligned_seq_qza = f'{self.workdir}/mafft-aligned-sequences-masked.qza'
        self.unrooted_tree_qza = f'{self.workdir}/unrooted-fasttree.qza'
        self.rooted_tree_qza = f'{self.workdir}/rooted-fasttree.qza'

        self.output_paths = (
            self.aligned_seq_qza,
            self.masked_aligned_seq_qza,
            self.unrooted_tree_qza,
            self.rooted_tree_qza
        )

    def execute(self):
        lines = [
            'qiime phylogeny align-to-tree-mafft-fasttree',
            f'--i-sequences {self.seq_qza}',
            f'--p-n-threads {self.threads}',
            f'--o-alignment {self.aligned_seq_qza}',
            f'--o-masked-alignment {self.masked_aligned_seq_qza}',
            f'--o-tree {self.unrooted_tree_qza}',
            f'--o-rooted-tree {self.rooted_tree_qza}'
        ]
        self.call(' \\\n  '.join(lines))
