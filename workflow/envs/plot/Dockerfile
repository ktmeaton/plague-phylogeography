FROM continuumio/miniconda3

LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for ktmeaton/plague-phylogeography: plot."

COPY plot.yaml .
RUN conda install -c conda-forge mamba \
    && conda env create -f plot.yaml \
    && conda clean -a \
    && echo "source activate plague-phylogeography-plot" > ~/.bashrc

ENV PATH /opt/conda/envs/plague-phylogeography-plot/bin:$PATH
