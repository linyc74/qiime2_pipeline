import os
from ete3 import Tree, TreeStyle, TextFace
from typing import Tuple
from .template import Processor
from .tools import edit_fpath
from .exporting import ExportAlignedSequence, ExportTree


class Phylogeny(Processor):

    seq_qza: str

    rooted_tree_qza: str
    rooted_tree_nwk: str
    unrooted_tree_nwk: str

    def main(self, seq_qza: str) -> str:
        self.seq_qza = seq_qza

        self.rooted_tree_qza, self.rooted_tree_nwk, self.unrooted_tree_nwk = \
            MafftFasttree(self.settings).main(seq_qza=self.seq_qza)

        DrawTree(self.settings).main(nwk=self.rooted_tree_nwk)
        DrawTree(self.settings).main(nwk=self.unrooted_tree_nwk)

        return self.rooted_tree_qza


class MafftFasttree(Processor):

    DSTDIR_NAME = 'phylogeny'

    seq_qza: str

    dstdir: str

    aligned_seq_qza: str
    masked_aligned_seq_qza: str
    unrooted_tree_qza: str
    rooted_tree_qza: str

    aligned_seq_fa: str
    masked_aligned_seq_fa: str
    unrooted_tree_nwk: str
    rooted_tree_nwk: str

    def main(self, seq_qza: str) -> Tuple[str, str, str]:
        self.seq_qza = seq_qza

        self.make_dstdir()
        self.set_output_paths()
        self.execute()
        self.export_to_fa()
        self.export_to_nwk()

        return self.rooted_tree_qza, self.rooted_tree_nwk, self.unrooted_tree_nwk

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def set_output_paths(self):
        self.aligned_seq_qza = f'{self.workdir}/mafft-aligned-sequences.qza'
        self.aligned_seq_fa = f'{self.dstdir}/mafft-aligned-sequences.fa'

        self.masked_aligned_seq_qza = f'{self.workdir}/mafft-aligned-sequences-masked.qza'
        self.masked_aligned_seq_fa = f'{self.dstdir}/mafft-aligned-sequences-masked.fa'

        self.rooted_tree_qza = f'{self.workdir}/fasttree-rooted.qza'
        self.rooted_tree_nwk = f'{self.dstdir}/fasttree-rooted.nwk'

        self.unrooted_tree_qza = f'{self.workdir}/fasttree-unrooted.qza'
        self.unrooted_tree_nwk = f'{self.dstdir}/fasttree-unrooted.nwk'

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

    def export_to_fa(self):
        for qza, fa in [
            (self.aligned_seq_qza, self.aligned_seq_fa),
            (self.masked_aligned_seq_qza, self.masked_aligned_seq_fa),
        ]:
            f = ExportAlignedSequence(self.settings).main(aligned_sequence_qza=qza)
            self.call(f'mv {f} {fa}')

    def export_to_nwk(self):
        for qza, nwk in [
            (self.rooted_tree_qza, self.rooted_tree_nwk),
            (self.unrooted_tree_qza, self.unrooted_tree_nwk),
        ]:
            f = ExportTree(self.settings).main(tree_qza=qza)
            self.call(f'mv {f} {nwk}')


class DrawTree(Processor):

    NEWICK_FORMAT = 1  # flexible with internal node names
    QUOTED_NODE_NAMES = True

    WIDTH = 200  # mm
    HEIGHT = 200
    DPI = 600

    def main(self, nwk: str):
        tree = Tree(
            nwk,
            format=self.NEWICK_FORMAT,
            quoted_node_names=self.QUOTED_NODE_NAMES)

        png = edit_fpath(
            fpath=nwk,
            old_suffix='.nwk',
            new_suffix='.png'
        )

        tree.render(
            file_name=png,
            w=self.WIDTH,
            h=self.HEIGHT,
            dpi=self.DPI,
            units='mm',
            tree_style=MyTreeStyle()
        )


class MyTreeStyle(TreeStyle):

    def __init__(self):
        super().__init__()
        self.show_leaf_name = False
        self.mode = 'c'
        self.arc_start = -180  # 0 degrees = 3 o'clock
        self.arc_span = 360
        self.title.add_face(TextFace("Hello ETE", fsize=20), column=0)
