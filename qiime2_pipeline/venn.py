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
        if n == 2 or n == 3:
            return True
        else:
            if n == 1:
                msg = f'There is only {n} group keyword: {self.group_keywords}, skip Venn diagram'
            else:  # n >= 4
                msg = f'There are {n} group keywords: {self.group_keywords}, skip Venn diagram'
            self.logger.info(msg)
            return False

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def plot_venn_diagrams(self):
        for tsv in self.tsvs:
            PlotOneVenn(self.settings).main(
                tsv=tsv,
                group_keywords=self.group_keywords,
                dstdir=self.dstdir)


class PlotOneVenn(Processor):

    FIGSIZE = (8, 6)
    DPI = 300

    tsv: str
    group_keywords: List[str]
    dstdir: str

    df: pd.DataFrame
    group_keyword_to_features: Dict[str, Set[str]]

    def main(
            self,
            tsv: str,
            group_keywords: List[str],
            dstdir: str):

        self.tsv = tsv
        self.group_keywords = group_keywords
        self.dstdir = dstdir

        self.assert_number_of_groups()
        self.read_tsv()
        self.init_group_keyword_to_features()
        self.count_features_for_each_group()

        self.init_figure()
        self.plot_venn()
        self.save_figure()

    def assert_number_of_groups(self):
        assert len(self.group_keywords) in [2, 3]

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

    def init_figure(self):
        plt.figure(figsize=self.FIGSIZE, dpi=self.DPI)

    def plot_venn(self):
        set_labels, subsets = [], []
        for group, features in self.group_keyword_to_features.items():
            set_labels.append(group)
            subsets.append(features)

        if len(subsets) == 2:
            venn2(subsets=subsets, set_labels=set_labels)
        else:
            venn3(subsets=subsets, set_labels=set_labels)

    def save_figure(self):
        png = edit_fpath(
            fpath=self.tsv,
            old_suffix='.tsv',
            new_suffix='.png',
            dstdir=self.dstdir)
        plt.savefig(png)
