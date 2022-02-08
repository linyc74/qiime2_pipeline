import argparse
import qiime2_pipeline


__VERSION__ = '1.6.1-beta'


PROG = 'python qiime2_pipeline'
DESCRIPTION = f'Custom-built Qiime2 pipeline (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'
REQUIRED = [
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
    {
        'keys': ['-2', '--fq2-suffix'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'suffix of read 2 fastq files',
        }
    },
    {
        'keys': ['-b', '--nb-classifier-qza'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'pre-trained naive Bayes classifier (.qza file) for feature classification (https://docs.qiime2.org/2021.8/data-resources)',
        }
    },
]
OPTIONAL = [
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
        'keys': ['-m', '--paired-end-mode'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'concat',
            'help': 'mode to combine paired end reads, "concat", "merge", or "pool" (default: %(default)s)',
        }
    },
    {
        'keys': ['-k', '--group-keywords'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'comma-separated group keywords, e.g. "control,treatment" (default: %(default)s)',
        }
    },
    {
        'keys': ['-i', '--otu-identity'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.97,
            'help': 'sequence identity (range 0, 1) for de novo OTU clustering (default: %(default)s)',
        }
    },
    {
        'keys': ['--classifier-reads-per-batch'],
        'properties': {
            'type': int,
            'required': False,
            'default': 0,
            'help': 'number of reads per batch for feature classifier, default indicates \'auto\' (default: %(default)s)',
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
        'keys': ['--skip-otu'],
        'properties': {
            'action': 'store_true',
            'help': 'use ASV without OTU clustering',
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
        qiime2_pipeline.Main().main(
            fq_dir=args.fq_dir,
            fq1_suffix=args.fq1_suffix,
            fq2_suffix=args.fq2_suffix,
            nb_classifier_qza=args.nb_classifier_qza,
            paired_end_mode=args.paired_end_mode,
            group_keywords=args.group_keywords,
            otu_identity=args.otu_identity,
            skip_otu=args.skip_otu,
            classifier_reads_per_batch=args.classifier_reads_per_batch,
            alpha_metrics=args.alpha_metrics,
            outdir=args.outdir,
            threads=args.threads,
            debug=args.debug)


if __name__ == '__main__':
    EntryPoint().main()
