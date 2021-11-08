FROM continuumio/anaconda3:2021.05

RUN wget https://data.qiime2.org/distro/core/qiime2-2021.8-py38-linux-conda.yml \
 && conda env create -n qiime2-2021.8 --file qiime2-2021.8-py38-linux-conda.yml \
 && rm qiime2-2021.8-py38-linux-conda.yml \
 && echo "conda activate qiime2-2021.8" >> /root/.bashrc

COPY ./qiime2_pipeline/* /qiime2_pipeline/qiime2_pipeline/
COPY ./__main__.py /qiime2_pipeline/
WORKDIR /
