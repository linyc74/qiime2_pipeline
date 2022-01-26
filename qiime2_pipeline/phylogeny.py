from .template import Processor, Settings
from .exporting import ExportAlignedSequence, ExportTree


class MafftFasttree(Processor):

    seq_qza: str

    aligned_seq_qza: str
    masked_aligned_seq_qza: str
    unrooted_tree_qza: str
    rooted_tree_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.export_aligned_sequence = ExportAlignedSequence(self.settings).main
        self.export_tree = ExportTree(self.settings).main

    def main(self, seq_qza: str) -> str:
        self.seq_qza = seq_qza

        self.set_output_paths()
        self.execute()
        self.export()

        return self.rooted_tree_qza

    def set_output_paths(self):
        self.aligned_seq_qza = f'{self.workdir}/mafft-aligned-sequences.qza'
        self.masked_aligned_seq_qza = f'{self.workdir}/mafft-aligned-sequences-masked.qza'
        self.unrooted_tree_qza = f'{self.workdir}/unrooted-fasttree.qza'
        self.rooted_tree_qza = f'{self.workdir}/rooted-fasttree.qza'

    def execute(self):
        cmd = self.CMD_LINEBREAK.join([
            'qiime phylogeny align-to-tree-mafft-fasttree',
            f'--i-sequences {self.seq_qza}',
            f'--p-n-threads {self.threads}',
            f'--o-alignment {self.aligned_seq_qza}',
            f'--o-masked-alignment {self.masked_aligned_seq_qza}',
            f'--o-tree {self.unrooted_tree_qza}',
            f'--o-rooted-tree {self.rooted_tree_qza}'
        ])
        self.call(cmd)

    def export(self):
        for aligned_sequence_qza in [self.aligned_seq_qza, self.masked_aligned_seq_qza]:
            fa = self.export_aligned_sequence(aligned_sequence_qza=aligned_sequence_qza)
            self.call(f'mv {fa} {self.outdir}/')

        for tree_qza in [self.rooted_tree_qza, self.unrooted_tree_qza]:
            nwk = self.export_tree(tree_qza=tree_qza)
            self.call(f'mv {nwk} {self.outdir}/')
