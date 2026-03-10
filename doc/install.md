# Install

```bash
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-linux-conda.yml
conda env create -n qiime2-amplicon-2024.2 --file qiime2-amplicon-2024.2-py38-linux-conda.yml
rm qiime2-amplicon-2024.2-py38-linux-conda.yml

conda activate qiime2-amplicon-2024.2

conda install -c bioconda trim-galore=0.6.6

pip install --no-cache-dir \
    lefse==1.1.2 \
    rpy2==3.5.10 \
    matplotlib-venn==0.11.7 \
    PyQt5==5.15.6 \
    ete3==3.1.2 \
    scikit-bio==0.5.6 \
    openpyxl==3.1.2 \
    venny4py==1.0.3

Rscript -e 'install.packages("survival", version="2.44", repos="https://cran.csie.ntu.edu.tw/")'
Rscript -e 'install.packages("mvtnorm", version="1.1", repos="https://cran.csie.ntu.edu.tw/")'
Rscript -e 'install.packages("modeltools", version="0.2", repos="https://cran.csie.ntu.edu.tw/")'
Rscript -e 'install.packages("coin", version="1.4", repos="https://cran.csie.ntu.edu.tw/")'
```

To activate environment and variables

```bash
conda activate qiime2-amplicon-2024.2 && export QT_QPA_PLATFORM=offscreen && export UNIFRAC_USE_GPU=N
```
