FROM continuumio/anaconda3:2021.11

RUN wget https://data.qiime2.org/distro/core/qiime2-2021.11-py38-linux-conda.yml \
 && conda env create -n qiime2-2021.11 --file qiime2-2021.11-py38-linux-conda.yml \
 && rm qiime2-2021.11-py38-linux-conda.yml

RUN conda install -c bioconda -n qiime2-2021.11 trim-galore=0.6.6

RUN conda install -c anaconda -n qiime2-2021.11 scikit-bio

ENV PATH /opt/conda/envs/qiime2-2021.11/bin:$PATH

COPY ./qiime2_pipeline/* /qiime2_pipeline/qiime2_pipeline/
COPY ./__main__.py /qiime2_pipeline/
WORKDIR /
