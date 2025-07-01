import argparse
import qiime2_pipeline


__VERSION__ = '2.9.6'


PROG = 'python qiime2_pipeline'
DESCRIPTION = f'Custom-built Qiime2 pipeline (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'
REQUIRED = [
    {
        'keys': ['-s', '--sample-sheet'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the sample sheet (CSV, TSV, or XLSX format), "Sample" and "Group" cloumns are required',
        }
    },
    {
        'keys': ['-f', '--fq-dir'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the directory containing all input fastq files',
        }
    },
    {
        'keys': ['-1', '--fq1-suffix'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'suffix of read 1 fastq files',
        }
    },
]
OPTIONAL = [
    {
        'keys': ['-2', '--fq2-suffix'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'suffix of read 2 fastq files, "None" for single end reads (default: %(default)s)',
        }
    },
    {
        'keys': ['-o', '--outdir'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'qiime2_pipeline_outdir',
            'help': 'path to the output directory (default: %(default)s)',
        }
    },
    {
        'keys': ['--sequencing-platform'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'illumina',
            'choices': ['illumina', 'pacbio', 'nanopore'],
            'help': 'sequencing platform (default: %(default)s)',
        }
    },
    {
        'keys': ['--clip-r1-5-prime'],
        'properties': {
            'type': int,
            'required': False,
            'default': 0,
            'help': 'hard clip <int> bp from 5\' end of Illumina read 1 (default: %(default)s)',
        }
    },
    {
        'keys': ['--clip-r2-5-prime'],
        'properties': {
            'type': int,
            'required': False,
            'default': 0,
            'help': 'hard clip <int> bp from 5\' end of Illumina read 2 (default: %(default)s)',
        }
    },
    {
        'keys': ['--paired-end-mode'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'merge',
            'choices': ['merge', 'pool'],
            'help': 'mode to combine Illumina paired end reads (default: %(default)s)',
        }
    },
    {
        'keys': ['--max-expected-error-bases'],
        'properties': {
            'type': float,
            'required': False,
            'default': 2.0,
            'help': 'max number of expected error bases for DADA2, i.e. the sum error rates across all bases (default: %(default)s)',
        }
    },
    {
        'keys': ['--otu-identity'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.97,
            'help': 'sequence identity (range 0, 1) for de novo OTU clustering (default: %(default)s)',
        }
    },
    {
        'keys': ['--skip-otu'],
        'properties': {
            'action': 'store_true',
            'help': 'use ASV without OTU clustering',
        }
    },
    {
        'keys': ['--feature-classifier'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'nb',
            'choices': ['nb', 'vsearch'],
            'help': 'taxonomy classification algorithm, "nb" for Naive Bayes (default: %(default)s)',
        }
    },
    {
        'keys': ['--nb-classifier-qza'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'pre-trained Naive Bayes classifier (.qza file) required for "nb" feature-classifier (default: %(default)s)',
        }
    },
    {
        'keys': ['--classifier-reads-per-batch'],
        'properties': {
            'type': int,
            'required': False,
            'default': 0,
            'help': 'reads per batch for "nb" feature-classifier, default indicates "auto" (default: %(default)s)',
        }
    },
    {
        'keys': ['--reference-sequence-qza'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'reference sequence (.qza file) required for "vsearch" feature-classifier (default: %(default)s)',
        }
    },
    {
        'keys': ['--reference-taxonomy-qza'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'reference taxonomy (.qza file) required for "vsearch" feature-classifier (default: %(default)s)',
        }
    },
    {
        'keys': ['--vsearch-classifier-max-hits'],
        'properties': {
            'type': str,
            'required': False,
            'default': 10,
            'help': 'maximum number of hits for the consensus of "vsearch" feature-classifier (default: %(default)s)',
        }
    },
    {
        'keys': ['--alpha-metrics'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'all',
            'help': 'comma-separated alpha-diversity metrics, default for all: '
                    '"chao1,shannon,gini_index,mcintosh_e,pielou_e,simpson,observed_features,fisher_alpha" (default: %(default)s)',
        }
    },
    {
        'keys': ['--beta-diversity-feature-level'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'feature',
            'choices': ['feature', 'species', 'genus', 'family', 'order', 'class', 'phylum'],
            'help': 'the level of features to be used for beta diversity (default: %(default)s)',
        }
    },
    {
        'keys': ['--heatmap-read-fraction'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.95,
            'help': 'fraction of total reads to be included in the heatmap (default: %(default)s)',
        }
    },
    {
        'keys': ['--n-taxa-barplot'],
        'properties': {
            'type': int,
            'required': False,
            'default': 20,
            'help': 'number of taxa shown in the percentage barplot (default: %(default)s)',
        }
    },
    {
        'keys': ['--colormap'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'Set1',
            'help': 'matplotlib colormap for plotting, or comma-separated color names, e.g. "darkred,lightgreen,skyblue" (default: %(default)s)',
        }
    },
    {
        'keys': ['--invert-colors'],
        'properties': {
            'action': 'store_true',
            'help': 'invert the order of colors',
        }
    },
    {
        'keys': ['--publication-figure'],
        'properties': {
            'action': 'store_true',
            'help': 'plot figures in the form and quality for paper publication',
        }
    },
    {
        'keys': ['--skip-differential-abundance'],
        'properties': {
            'action': 'store_true',
            'help': 'skip differential abundance analysis',
        }
    },
    {
        'keys': ['--differential-abundance-p-value'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.05,
            'help': 'p value cutoff for differential abundance (default: %(default)s)',
        }
    },
    {
        'keys': ['--run-picrust2'],
        'properties': {
            'action': 'store_true',
            'help': 'run PICRUSt2 analysis',
        }
    },
    {
        'keys': ['-t', '--threads'],
        'properties': {
            'type': int,
            'required': False,
            'default': 4,
            'help': 'number of CPU threads (default: %(default)s)',
        }
    },
    {
        'keys': ['-d', '--debug'],
        'properties': {
            'action': 'store_true',
            'help': 'debug mode',
        }
    },
    {
        'keys': ['-h', '--help'],
        'properties': {
            'action': 'help',
            'help': 'show this help message',
        }
    },
    {
        'keys': ['-v', '--version'],
        'properties': {
            'action': 'version',
            'version': __VERSION__,
            'help': 'show version',
        }
    },
]


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        self.run()

    def set_parser(self):
        self.parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        for item in OPTIONAL:
            group.add_argument(*item['keys'], **item['properties'])

    def run(self):
        args = self.parser.parse_args()
        print(f'Start running Qiime2 pipeline version {__VERSION__}\n', flush=True)
        qiime2_pipeline.main(
            sample_sheet=args.sample_sheet,
            fq_dir=args.fq_dir,
            fq1_suffix=args.fq1_suffix,
            fq2_suffix=args.fq2_suffix,
            outdir=args.outdir,

            sequencing_platform=args.sequencing_platform,

            clip_r1_5_prime=args.clip_r1_5_prime,
            clip_r2_5_prime=args.clip_r2_5_prime,

            paired_end_mode=args.paired_end_mode,
            max_expected_error_bases=args.max_expected_error_bases,

            otu_identity=args.otu_identity,
            skip_otu=args.skip_otu,

            feature_classifier=args.feature_classifier,
            nb_classifier_qza=args.nb_classifier_qza,
            classifier_reads_per_batch=args.classifier_reads_per_batch,
            reference_sequence_qza=args.reference_sequence_qza,
            reference_taxonomy_qza=args.reference_taxonomy_qza,
            vsearch_classifier_max_hits=args.vsearch_classifier_max_hits,

            alpha_metrics=args.alpha_metrics,
            beta_diversity_feature_level=args.beta_diversity_feature_level,
            heatmap_read_fraction=args.heatmap_read_fraction,
            n_taxa_barplot=args.n_taxa_barplot,
            colormap=args.colormap,
            invert_colors=args.invert_colors,
            publication_figure=args.publication_figure,
            skip_differential_abundance=args.skip_differential_abundance,
            differential_abundance_p_value=args.differential_abundance_p_value,
            run_picrust2=args.run_picrust2,

            threads=args.threads,
            debug=args.debug)


if __name__ == '__main__':
    EntryPoint().main()
