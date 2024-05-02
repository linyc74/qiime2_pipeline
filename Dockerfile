FROM continuumio/miniconda3:23.5.2-0

RUN apt-get update \
 && apt-get install -y \
    libgl1-mesa-glx \
    libx11-xcb1 \
    libfontconfig \
    libxkbcommon0 \
    libxkbcommon-x11-0 \
    libdbus-1-3 \
    build-essential

ARG yml="qiime2-2023.2-py38-linux-conda.yml"

RUN wget https://data.qiime2.org/distro/core/${yml} \
 && conda env create -n qiime2 --file ${yml} \
 && rm ${yml}

ARG bin_path="/opt/conda/envs/qiime2/bin"

ENV PATH ${bin_path}:$PATH

# rpy2==3.5.10 to avoid lefse bug
RUN ${bin_path}/pip install --no-cache-dir \
    lefse==1.1.2 \
    rpy2==3.5.10 \
    matplotlib-venn==0.11.7 \
    PyQt5==5.15.6 \
    ete3==3.1.2 \
    scikit-bio==0.5.6

ENV QT_QPA_PLATFORM offscreen

RUN conda install -c bioconda -n qiime2 \
    trim-galore=0.6.6 \
 && conda install -c conda-forge \
    zip=3.0

RUN Rscript -e 'install.packages("survival", version="2.44", repos="https://cran.csie.ntu.edu.tw/")' \
 && Rscript -e 'install.packages("mvtnorm", version="1.1", repos="https://cran.csie.ntu.edu.tw/")' \
 && Rscript -e 'install.packages("modeltools", version="0.2", repos="https://cran.csie.ntu.edu.tw/")' \
 && Rscript -e 'install.packages("coin", version="1.4", repos="https://cran.csie.ntu.edu.tw/")'

RUN conda install -c bioconda \
    epa-ng=0.3.8 \
    gappa=0.7.1 \
 && Rscript -e 'install.packages("castor", version="1.7.11", repos="https://cran.csie.ntu.edu.tw/")' \
 && wget https://github.com/picrust/picrust2/archive/v2.5.2.tar.gz \
 && tar xvzf v2.5.2.tar.gz \
 && rm v2.5.2.tar.gz \
 && cd picrust2-2.5.2 \
 && pip install --editable . \
 && cd ..

RUN ${bin_path}/pip install --no-cache-dir openpyxl==3.1.2

RUN conda clean --all --yes

COPY ./qiime2_pipeline/* /qiime2_pipeline/qiime2_pipeline/
COPY ./__main__.py /qiime2_pipeline/
WORKDIR /
