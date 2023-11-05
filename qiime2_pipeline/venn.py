import os
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Dict, Set
from matplotlib_venn import venn2, venn3
from .tools import edit_fpath
from .template import Processor


GROUP_COLUMN = 'Group'


class PlotVennDiagrams(Processor):

    DSTDIR_NAME = 'venn'

    tsvs: List[str]
    sample_sheet: str

    dstdir: str

    def main(
            self,
            tsvs: List[str],
            sample_sheet: str):

        self.tsvs = tsvs
        self.sample_sheet = sample_sheet

        valid = self.check_number_of_groups()
        if not valid:
            return

        self.make_dstdir()
        self.plot_venn_diagrams()

    def check_number_of_groups(self) -> bool:
        groups = pd.read_csv(self.sample_sheet, index_col=0)[GROUP_COLUMN].unique()

        n = len(groups)
        if n in [2, 3]:
            return True
        else:
            p = 'There is only 1' if n == 1 else f'There are {n}'
            msg = f'{p} group keywords: {groups}, skip Venn diagram'
            self.logger.info(msg)
            return False

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def plot_venn_diagrams(self):
        for tsv in self.tsvs:
            ProcessTsvPlotVenn(self.settings).main(
                tsv=tsv,
                sample_sheet=self.sample_sheet,
                dstdir=self.dstdir)


class ProcessTsvPlotVenn(Processor):

    FIGSIZE = (8, 6)
    DPI = 300

    tsv: str
    sample_sheet: str
    dstdir: str

    df: pd.DataFrame
    sample_to_group: Dict[str, str]
    group_to_features: Dict[str, Set[str]]

    def main(
            self,
            tsv: str,
            sample_sheet: str,
            dstdir: str):

        self.tsv = tsv
        self.sample_sheet = sample_sheet
        self.dstdir = dstdir

        self.read_tsv()
        self.set_sample_to_group()
        self.init_group_to_features()
        self.count_features_for_each_group()
        self.plot_venn()

    def read_tsv(self):
        self.df = pd.read_csv(self.tsv, sep='\t', index_col=0)

    def set_sample_to_group(self):
        df = pd.read_csv(self.sample_sheet, index_col=0)
        self.sample_to_group = df[GROUP_COLUMN].to_dict()

    def init_group_to_features(self):
        df = pd.read_csv(self.sample_sheet, index_col=0)
        groups = df[GROUP_COLUMN].unique()
        self.group_to_features = {g: set() for g in groups}

    def count_features_for_each_group(self):
        for sample in self.df.columns:
            group = self.sample_to_group[sample]
            set_1 = self.group_to_features[group]

            set_2 = self.df[sample][self.df[sample] > 0].index

            self.group_to_features[group] = set_1.union(set_2)

    def plot_venn(self):
        groups = list(self.group_to_features.keys())
        subsets = [self.group_to_features[g] for g in groups]

        output_prefix = edit_fpath(
            fpath=self.tsv,
            old_suffix='.tsv',
            new_suffix='',
            dstdir=self.dstdir)

        PlotVenn(self.settings).main(
            set_labels=groups,
            subsets=subsets,
            output_prefix=output_prefix)


class PlotVenn(Processor):

    FIGSIZE = (8, 6)
    DPI = 600

    set_labels: List[str]
    subsets: List[Set[str]]
    output_prefix: str

    def main(
            self,
            set_labels: List[str],
            subsets: List[Set[str]],
            output_prefix: str):

        self.set_labels = set_labels
        self.subsets = subsets
        self.output_prefix = output_prefix

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
        for ext in ['pdf', 'png']:
            plt.savefig(f'{self.output_prefix}.{ext}', dpi=self.DPI)
        plt.close()
