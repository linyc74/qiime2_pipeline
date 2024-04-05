import os
from random import randrange
from typing import Tuple, Dict, List
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, RectFace
from .tools import edit_fpath
from .template import Processor
from .exporting import ExportAlignedSequence, ExportTree


class Phylogeny(Processor):

    DSTDIR_NAME = 'phylogeny'

    seq_qza: str

    dstdir: str

    rooted_tree_qza: str
    rooted_tree_nwk: str
    unrooted_tree_nwk: str

    def main(self, seq_qza: str) -> str:

        self.seq_qza = seq_qza

        self.make_dstdir()
        self.mafft_fasttree()
        self.draw_trees()

        return self.rooted_tree_qza

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def mafft_fasttree(self):
        self.rooted_tree_qza, self.rooted_tree_nwk, self.unrooted_tree_nwk = \
            MafftFasttree(self.settings).main(
                seq_qza=self.seq_qza,
                dstdir=self.dstdir)

    def draw_trees(self):
        for nwk in [self.rooted_tree_nwk, self.unrooted_tree_nwk]:
            DrawTree(self.settings).main(nwk=nwk, dstdir=self.dstdir)


class MafftFasttree(Processor):

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

    def main(self, seq_qza: str, dstdir: str) -> Tuple[str, str, str]:

        self.seq_qza = seq_qza
        self.dstdir = dstdir

        self.set_output_paths()
        self.execute()
        self.export_to_fa()
        self.export_to_nwk()

        return self.rooted_tree_qza, self.rooted_tree_nwk, self.unrooted_tree_nwk

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
        log = f'{self.outdir}/qiime-phylogeny-align-to-tree-mafft-fasttree.log'
        cmd = self.CMD_LINEBREAK.join([
            'qiime phylogeny align-to-tree-mafft-fasttree',
            f'--i-sequences {self.seq_qza}',
            f'--p-n-threads {self.threads}',
            f'--o-alignment {self.aligned_seq_qza}',
            f'--o-masked-alignment {self.masked_aligned_seq_qza}',
            f'--o-tree {self.unrooted_tree_qza}',
            f'--o-rooted-tree {self.rooted_tree_qza}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ])
        self.call(cmd)

    def export_to_fa(self):
        for qza, fa in [
            (self.aligned_seq_qza, self.aligned_seq_fa),
            (self.masked_aligned_seq_qza, self.masked_aligned_seq_fa),
        ]:
            f = ExportAlignedSequence(self.settings).main(aligned_sequence_qza=qza)
            self.call(f'mv "{f}" "{fa}"')

    def export_to_nwk(self):
        for qza, nwk in [
            (self.rooted_tree_qza, self.rooted_tree_nwk),
            (self.unrooted_tree_qza, self.unrooted_tree_nwk),
        ]:
            f = ExportTree(self.settings).main(tree_qza=qza)
            self.call(f'mv "{f}" "{nwk}"')


class DrawTree(Processor):

    NEWICK_FORMAT = 1  # flexible with internal node names
    QUOTED_NODE_NAMES = True
    FIG_WIDTH = 200  # mm
    DPI = 1200
    LINE_WIDTH = 1
    LEAF_LENGTH = 100
    COLORS = [
        'RoyalBlue',
        'Crimson',
        'DarkSlateBlue',
        'DarkOrange',
        'Navy',
        'MediumVioletRed',
        'Indigo',
        'DodgerBlue',
        'DarkMagenta',
        'DarkGreen',
        'SteelBlue',
        'DarkRed',
        'Teal',
        'MediumBlue',
        'MidnightBlue',
        'SaddleBrown',
        'Maroon',
        'DeepSkyBlue',
        'ForestGreen',
    ]

    nwk: str
    dstdir: str

    tree: Tree
    sorted_phylums: List[str]
    phylum_to_color: Dict[str, str]

    def main(self, nwk: str, dstdir: str):
        self.nwk = nwk
        self.dstdir = dstdir

        self.set_tree()
        self.sort_phylums_by_count()
        self.set_phylum_to_color()
        self.traverse_tree_to_color_leafs()
        self.render_tree()

    def set_tree(self):
        self.tree = Tree(
            self.nwk,
            format=self.NEWICK_FORMAT,
            quoted_node_names=self.QUOTED_NODE_NAMES)

    def sort_phylums_by_count(self):
        leaf_phylums = map(
            leaf_name_to_bacterial_phylum, self.tree.get_leaf_names()
        )
        phylum_to_count = {}
        for p in leaf_phylums:
            phylum_to_count.setdefault(p, 0)
            phylum_to_count[p] += 1

        def __get_count(tup: Tuple[str, int]) -> int:
            return tup[1]

        sorted_phylum_and_count = sorted(
            phylum_to_count.items(),
            key=__get_count,
            reverse=True
        )
        self.sorted_phylums, _ = zip(*sorted_phylum_and_count)  # unzip

    def set_phylum_to_color(self):
        self.phylum_to_color = {}
        for i, phylum in enumerate(self.sorted_phylums):
            self.phylum_to_color[phylum] = self.COLORS[i] if i < len(self.COLORS) else random_color()

    def traverse_tree_to_color_leafs(self):
        no_node = NodeStyle(
            size=0,
            vt_line_width=self.LINE_WIDTH,
            hz_line_width=self.LINE_WIDTH)

        for node in self.tree.traverse():
            node.set_style(no_node)
            if node.is_leaf():
                phylum = leaf_name_to_bacterial_phylum(node.name)
                color = self.phylum_to_color[phylum]
                face = RectFace(
                    width=self.LEAF_LENGTH,
                    height=self.LINE_WIDTH,
                    fgcolor=color,
                    bgcolor=color)
                node.add_face(face=face, column=0)

    def render_tree(self):
        for fmt in ['png', 'pdf', 'svg']:
            file_name = edit_fpath(
                fpath=self.nwk,
                old_suffix='.nwk',
                new_suffix=f'.{fmt}',
                dstdir=self.dstdir)

            tree_style = MyTreeStyle(
                phylum_to_color=self.phylum_to_color)

            self.tree.render(
                file_name=file_name,
                w=self.FIG_WIDTH,
                dpi=self.DPI,
                units='mm',
                tree_style=tree_style)


class MyTreeStyle(TreeStyle):

    def __init__(self, phylum_to_color: Dict[str, str]):
        super().__init__()

        self.show_leaf_name = False
        self.mode = 'c'
        self.arc_start = -180  # 0 degrees = 3 o'clock
        self.arc_span = 360
        self.root_opening_factor = 0
        self.min_leaf_separation = 0

        self.set_phylum_color_legends(phylum_to_color)

    def set_phylum_color_legends(self, phylum_to_color):
        for phylum, color in phylum_to_color.items():
            t = TextFace(text=phylum, fgcolor=color, ftype='Arial')
            self.legend.add_face(t, column=1)


def leaf_name_to_bacterial_phylum(leaf_name: str) -> str:
    """
    'ASV_0053; d__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales' -> 'Bacteroidetes'
    """
    if '__Bacteria;' not in leaf_name:
        return 'Unassigned'
    if '; p__' not in leaf_name:
        return 'Unassigned'
    else:
        return leaf_name.split('; p__')[1].split(';')[0]


def random_color() -> str:
    rgb = tuple(randrange(0, 150) for _ in range(3))
    return '#%02x%02x%02x' % rgb
