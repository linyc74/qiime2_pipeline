[![DOI](https://zenodo.org/badge/425668537.svg)](https://zenodo.org/badge/latestdoi/425668537)

# Qiime2 Pipeline

**Custom-built Qiime2 pipeline for 16S rRNA sequence analysis**

## Usage

```bash
git clone https://github.com/linyc74/qiime2_pipeline.git

python qiime2_pipeline \
  -f FQ_DIR \
  -1 FQ1_SUFFIX \
  -2 FQ2_SUFFIX \
  -b NB_CLASSIFIER_QZA
```

The program automatically detects all fastq files in the `FQ_DIR` directory,
using suffixes of read 1 (`FQ1_SUFFIX`) and read 2 (`FQ2_SUFFIX`) files,
then run the computational pipeline. Taxonomic classification is performed
using pre-trained naïve Bayes classifier (`NB_CLASSIFIER_QZA`), which can be
downloaded [here](https://data.qiime2.org/2022.8/common/silva-138-99-nb-classifier.qza)
from the [Qiime2 data resources page](https://docs.qiime2.org/2022.8/data-resources/).

For more options, see help message by

```bash
python qiime2_pipeline --help
```

## Environment

Assuming [Anaconda](https://www.anaconda.com/) has already been installed, create a Qiime2 environment named `qiime2`:

```bash
wget https://data.qiime2.org/distro/core/qiime2-2021.11-py38-linux-conda.yml
conda env create --name qiime2 --file qiime2-2021.11-py38-linux-conda.yml
```

Activate the `qiime2` environment:

```bash
conda activate qiime2
```

Additional packages to be installed:

```bash
conda install --channel bioconda --name qiime2 trim-galore

# Python packages
conda install --channel anaconda --name qiime2 scikit-bio
pip install lefse matplotlib-venn PyQt5 ete3

# R packages
Rscript -e 'install.packages("survival", version="2.44", repos="https://cran.csie.ntu.edu.tw/")'
Rscript -e 'install.packages("mvtnorm", version="1.1", repos="https://cran.csie.ntu.edu.tw/")'
Rscript -e 'install.packages("modeltools", version="0.2", repos="https://cran.csie.ntu.edu.tw/")'
Rscript -e 'install.packages("coin", version="1.4", repos="https://cran.csie.ntu.edu.tw/")'
```

## Docker

Pull the docker image, then run the python command in the container:

```bash
docker pull linyc74/qiime2-pipeline:latest

docker run \
  --volume "path/to/fq_dir":"path/to/fq_dir" \
  linyc74/qiime2-pipeline:latest \
  python qiime2_pipeline \
    -f FQ_DIR \
    -1 FQ1_SUFFIX \
    -2 FQ2_SUFFIX \
    -b NB_CLASSIFIER_QZA
```
