import argparse
import qiime2_pipeline


__version__ = '1.0.0-beta'


PROG = 'python qiime2_pipeline'
DESCRIPTION = f'Custom built Qiime2 pipeline (version {__version__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'
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
    {
        'keys': ['-l', '--read-length'],
        'properties': {
            'type': int,
            'required': True,
            'help': 'read length (bp)',
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
            'version': __version__,
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
        qiime2_pipeline.Main().main(
            fq_dir=args.fq_dir,
            fq1_suffix=args.fq1_suffix,
            fq2_suffix=args.fq2_suffix,
            nb_classifier_qza=args.nb_classifier_qza,
            read_length=args.read_length,
            outdir=args.outdir,
            threads=args.threads,
            debug=args.debug)


if __name__ == '__main__':
    EntryPoint().main()