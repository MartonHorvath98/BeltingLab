FROM condaforge/miniforge3:latest

LABEL image.author.name="MartonHorvath98"
LABEL image.author.email="marton.hrvth.98@gmail.com"

COPY ./env.yml /tmp/env.yml

RUN mamba env create -f /tmp/env.yml && \
    mamba clean --all --yes

ENV PATH=/opt/conda/envs/rnaseq-nf/bin:$PATH
