from os.path import basename
from .template import Processor, Settings


class MafftFasttree(Processor):

    seq_qza: str

    aligned_seq_qza: str
    masked_aligned_seq_qza: str
    unrooted_tree_qza: str
    rooted_tree_qza: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

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

    def export(self):
        for qza in [self.aligned_seq_qza, self.masked_aligned_seq_qza]:
            ExportAlignedSequence(self.settings).main(
                aligned_sequence_qza=qza)

        for tree_qza in [self.rooted_tree_qza, self.unrooted_tree_qza]:
            ExportTree(self.settings).main(
                tree_qza=tree_qza)


class ExportAlignedSequence(Processor):

    aligned_sequence_qza: str
    output_fa: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, aligned_sequence_qza: str):

        self.aligned_sequence_qza = aligned_sequence_qza

        self.qza_to_fa()
        self.set_output_fa()
        self.move_fa()

    def qza_to_fa(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.aligned_sequence_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def set_output_fa(self):
        fname = basename(self.aligned_sequence_qza)[:-len('.qza')] + '.fa'
        self.output_fa = f'{self.outdir}/{fname}'

    def move_fa(self):
        self.call(f'mv {self.workdir}/aligned-dna-sequences.fasta {self.output_fa}')


class ExportTree(Processor):

    tree_qza: str
    output_nwk: str

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def main(self, tree_qza: str):

        self.tree_qza = tree_qza

        self.qza_to_nwk()
        self.set_output_nwk(tree_qza)
        self.move_nwk()

    def qza_to_nwk(self):
        lines = [
            'qiime tools export',
            f'--input-path {self.tree_qza}',
            f'--output-path {self.workdir}',
        ]
        self.call(' \\\n  '.join(lines))

    def set_output_nwk(self, tree_qza):
        fname = basename(tree_qza)[:-len('.qza')] + '.nwk'
        self.output_nwk = f'{self.outdir}/{fname}'

    def move_nwk(self):
        self.call(f'mv {self.workdir}/tree.nwk {self.output_nwk}')
