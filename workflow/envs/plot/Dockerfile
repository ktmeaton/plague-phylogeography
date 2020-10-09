FROM continuumio/miniconda3
LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for ktmeaton/plague-phylogeography: plot."
COPY plot.yaml /

# install mamba and create environment
RUN conda install -c conda-forge mamba; \
    mamba env create -f /plot.yaml;

# Cleanup conda to reduce image size
RUN conda clean -a;

# Activate conda environment on login
RUN echo "source activate plague-phylogeography-plot" > ~/.bashrc;

ENV BASH_ENV ~/.bashrc
SHELL ["/bin/bash", "-c"]
