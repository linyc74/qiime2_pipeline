from .setup import TestCase
from qiime2_pipeline.beta import BetaDiversity
from qiime2_pipeline.denoise import Dada2PairedEnd, Dada2SingleEnd
from qiime2_pipeline.concat import Concat, BatchConcat, Pool, BatchPool
from qiime2_pipeline.importing import ImportPairedEndFastq, ImportSingleEndFastq


class MyTest(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def __test_concat(self):
        actual = Concat(self.settings).main(
            fq1=f'{self.indir}/R1.fastq.gz',
            fq2=f'{self.indir}/R2.fastq.gz'
        )
        expected = f'{self.workdir}/concat.fq'
        self.assertFileExists(expected, actual)

    def __test_batch_concat(self):
        fq_dir, fq_suffix = BatchConcat(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz')
        self.assertFileExists(f'{self.workdir}/concat_fastqs', fq_dir)
        self.assertEqual('.fq', fq_suffix)

    def __test_pool(self):
        actual = Pool(self.settings).main(
            fq1=f'{self.indir}/R1.fastq.gz',
            fq2=f'{self.indir}/R2.fastq.gz'
        )
        expected = f'{self.workdir}/pool.fq'
        self.assertFileExists(expected, actual)

    def __test_batch_pool(self):
        fq_dir, fq_suffix = BatchPool(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz')
        self.assertFileExists(f'{self.workdir}/pool_fastqs', fq_dir)
        self.assertEqual('.fq', fq_suffix)

    def __test_import_single_end_fastq(self):
        actual = ImportSingleEndFastq(self.settings).main(
            fq_dir=f'{self.indir}/concat_fastqs',
            fq_suffix='.fq')
        expected = f'{self.workdir}/single-end-demultiplexed.qza'
        self.assertFileExists(expected, actual)

    def __test_import_paired_end_fastq(self):
        actual = ImportPairedEndFastq(self.settings).main(
            fq_dir=f'{self.indir}/fq_dir',
            fq1_suffix='_L001_R1_001.fastq.gz',
            fq2_suffix='_L001_R2_001.fastq.gz')
        expected = f'{self.workdir}/paired-end-demultiplexed.qza'
        self.assertFileExists(expected, actual)

    def __test_dada2_single_end(self):
        Dada2SingleEnd(self.settings).main(
            demultiplexed_seq_qza=f'{self.indir}/single-end-demultiplexed.qza'
        )

    def __test_dada2_paired_end(self):
        Dada2PairedEnd(self.settings).main(
            demultiplexed_seq_qza=f'{self.indir}/paired-end-demultiplexed.qza'
        )

    def __test_beta_diversity(self):
        actual = BetaDiversity(self.settings).main(
            feature_table_qza=f'{self.indir}/feature-table.qza',
            rooted_tree_qza=f'{self.indir}/rooted-tree.qza'
        )
        expected = [
            f'{self.outdir}/beta-diversity/jaccard.tsv',
            f'{self.outdir}/beta-diversity/euclidean.tsv',
            f'{self.outdir}/beta-diversity/cosine.tsv',
            f'{self.outdir}/beta-diversity/weighted_normalized_unifrac.tsv',
            f'{self.outdir}/beta-diversity/weighted_unifrac.tsv',
            f'{self.outdir}/beta-diversity/generalized_unifrac.tsv',
            f'{self.outdir}/beta-diversity/unweighted_unifrac.tsv',
        ]
        for e, a in zip(expected, actual):
            self.assertFileExists(e, a)
