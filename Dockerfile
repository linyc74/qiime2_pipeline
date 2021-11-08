FROM continuumio/anaconda3:2021.05

COPY ./qiime2_pipeline/* /qiime2_pipeline/qiime2_pipeline/
COPY ./__main__.py /qiime2_pipeline/
WORKDIR /
