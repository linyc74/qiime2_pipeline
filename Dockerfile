FROM continuumio/anaconda3:2021.11

RUN wget https://data.qiime2.org/distro/core/qiime2-2021.11-py38-linux-conda.yml \
 && conda env create -n qiime2-2021.11 --file qiime2-2021.11-py38-linux-conda.yml \
 && rm qiime2-2021.11-py38-linux-conda.yml

RUN conda install -c bioconda -n qiime2-2021.11 trim-galore=0.6.6

RUN conda install -c anaconda -n qiime2-2021.11 scikit-bio=0.5.6

ENV PATH /opt/conda/envs/qiime2-2021.11/bin:$PATH

RUN /opt/conda/envs/qiime2-2021.11/bin/pip install --no-cache-dir \
    lefse==1.1.2 \
    matplotlib-venn==0.11.7

RUN Rscript -e 'install.packages("survival", version="2.44", repos="https://cran.csie.ntu.edu.tw/")' \
 && Rscript -e 'install.packages("mvtnorm", version="1.1", repos="https://cran.csie.ntu.edu.tw/")' \
 && Rscript -e 'install.packages("modeltools", version="0.2", repos="https://cran.csie.ntu.edu.tw/")' \
 && Rscript -e 'install.packages("coin", version="1.4", repos="https://cran.csie.ntu.edu.tw/")'

RUN apt-get update \
 && apt-get install -y \
    libgl1-mesa-glx \
    libx11-xcb1 \
    libfontconfig \
    libxkbcommon0 \
    libxkbcommon-x11-0 \
    libdbus-1-3 \
 && /opt/conda/envs/qiime2-2021.11/bin/pip install --no-cache-dir \
    PyQt5==5.15.6 \
    ete3==3.1.2
ENV QT_QPA_PLATFORM offscreen

RUN apt-get install -y --fix-missing zip

COPY ./qiime2_pipeline/* /qiime2_pipeline/qiime2_pipeline/
COPY ./__main__.py /qiime2_pipeline/
WORKDIR /
