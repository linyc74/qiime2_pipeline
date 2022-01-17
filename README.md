# Qiime2 Pipeline

**Custom-built Qiime2 pipeline**

## Usage

```bash
git clone https://github.com/linyc74/qiime2_pipeline.git

python qiime2_pipeline \
  -f FQ_DIR \
  -1 FQ1_SUFFIX \
  -2 FQ2_SUFFIX \
  -b NB_CLASSIFIER_QZA
```

For more options, see help message by

```bash
python qiime2_pipeline --help
```

## Environment

Create a Qiime2 environment:

```bash
wget https://data.qiime2.org/distro/core/qiime2-2021.11-py38-linux-conda.yml
conda env create -n qiime2-2021.11 --file qiime2-2021.11-py38-linux-conda.yml
rm qiime2-2021.11-py38-linux-conda.yml
```

Install TrimGalore:

```bash
conda install -c bioconda -n qiime2-2021.11 trim-galore
```

Activate the environment before usage.

```bash
conda activate qiime2-2021.11
```
