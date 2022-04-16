import os
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Dict, Set
from matplotlib_venn import venn2, venn3
from .tools import edit_fpath
from .template import Processor


class PlotVennDiagrams(Processor):

    DSTDIR_NAME = 'venn'

    tsvs: List[str]
    group_keywords: List[str]

    dstdir: str

    def main(
            self,
            tsvs: List[str],
            group_keywords: List[str]):

        self.tsvs = tsvs
        self.group_keywords = group_keywords

        valid = self.check_number_of_groups()
        if not valid:
            return

        self.make_dstdir()
        self.plot_venn_diagrams()

    def check_number_of_groups(self) -> bool:
        n = len(self.group_keywords)
        if n in [2, 3]:
            return True
        else:
            p = 'There is only 1' if n == 1 else f'There are {n}'
            msg = f'{p} group keywords: {self.group_keywords}, skip Venn diagram'
            self.logger.info(msg)
            return False

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def plot_venn_diagrams(self):
        for tsv in self.tsvs:
            ProcessTsvPlotVenn(self.settings).main(
                tsv=tsv,
                group_keywords=self.group_keywords,
                dstdir=self.dstdir)


class ProcessTsvPlotVenn(Processor):

    FIGSIZE = (8, 6)
    DPI = 300

    tsv: str
    group_keywords: List[str]
    dstdir: str

    df: pd.DataFrame
    group_keyword_to_features: Dict[str, Set[str]]
    png: str

    def main(
            self,
            tsv: str,
            group_keywords: List[str],
            dstdir: str):

        self.tsv = tsv
        self.group_keywords = group_keywords
        self.dstdir = dstdir

        self.read_tsv()
        self.init_group_keyword_to_features()
        self.count_features_for_each_group()
        self.set_png()
        self.plot_venn()

    def read_tsv(self):
        self.df = pd.read_csv(self.tsv, sep='\t', index_col=0)

    def init_group_keyword_to_features(self):
        self.group_keyword_to_features = {k: set() for k in self.group_keywords}

    def count_features_for_each_group(self):
        for c in self.df.columns:
            for k in self.group_keyword_to_features.keys():
                if k in c:
                    set_1 = self.group_keyword_to_features[k]
                    set_2 = self.df[c][self.df[c] > 0].index
                    self.group_keyword_to_features[k] = set_1.union(set_2)

    def set_png(self):
        self.png = edit_fpath(
            fpath=self.tsv,
            old_suffix='.tsv',
            new_suffix='.png',
            dstdir=self.dstdir)

    def plot_venn(self):
        subsets = [
            self.group_keyword_to_features[k] for k in self.group_keywords
        ]
        PlotVenn(self.settings).main(
            set_labels=self.group_keywords,
            subsets=subsets,
            png=self.png)


class PlotVenn(Processor):

    FIGSIZE = (8, 6)
    DPI = 300

    set_labels: List[str]
    subsets: List[Set[str]]
    png: str

    def main(
            self,
            set_labels: List[str],
            subsets: List[Set[str]],
            png: str):

        self.set_labels = set_labels
        self.subsets = subsets
        self.png = png

        self.assert_number_of_groups()
        self.init_figure()
        self.plot()
        self.save_figure()

    def assert_number_of_groups(self):
        assert len(self.set_labels) in [2, 3]

    def init_figure(self):
        plt.figure(figsize=self.FIGSIZE, dpi=self.DPI)

    def plot(self):
        if len(self.subsets) == 2:
            venn2(subsets=self.subsets, set_labels=self.set_labels)
        else:
            venn3(subsets=self.subsets, set_labels=self.set_labels)

    def save_figure(self):
        plt.savefig(self.png)
