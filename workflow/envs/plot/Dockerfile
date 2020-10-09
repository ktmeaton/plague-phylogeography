FROM continuumio/miniconda3
LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for plague-phylogeography: plot."
COPY plot.yaml /
RUN conda install -c conda-forge mamba && mamba env create -f /plot.yaml && conda clean -a
ENV PATH /opt/conda/envs/plague-phylogeography-plot/bin:$PATH
