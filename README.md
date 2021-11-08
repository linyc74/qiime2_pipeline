# Qiime2 Pipeline

**Custom-built Qiime2 pipeline**

```
python qiime2_pipeline \
  -f FQ_DIR \
  -1 FQ1_SUFFIX \
  -2 FQ2_SUFFIX \
  -b NB_CLASSIFIER_QZA \
  -l READ_LENGTH \
  -o OUTDIR \
  -t THREADS

required arguments:
  -f FQ_DIR, --fq-dir FQ_DIR
                        path to the directory containing all input fastq files
  -1 FQ1_SUFFIX, --fq1-suffix FQ1_SUFFIX
                        suffix of read 1 fastq files
  -2 FQ2_SUFFIX, --fq2-suffix FQ2_SUFFIX
                        suffix of read 2 fastq files
  -b NB_CLASSIFIER_QZA, --nb-classifier-qza NB_CLASSIFIER_QZA
                        pre-trained naive Bayes classifier (.qza file) for feature classification (https://docs.qiime2.org/2021.8/data-resources)
  -l READ_LENGTH, --read-length READ_LENGTH
                        read length (bp)

optional arguments:
  -o OUTDIR, --outdir OUTDIR
                        path to the output directory (default: qiime2_pipeline_outdir)
  -t THREADS, --threads THREADS
                        number of CPU threads (default: 4)
  -d, --debug           debug mode
  -h, --help            show this help message
  -v, --version         show version
```
