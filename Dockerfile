FROM continuumio/miniconda3

COPY ./ ./

RUN conda env create -f conda/env.yml
