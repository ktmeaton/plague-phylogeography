FROM continuumio/miniconda3

LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for ktmeaton/plague-phylogeography: qc."

COPY qc.yaml .
RUN conda install -c conda-forge mamba \
    && conda env create -f qc.yaml \
    && conda clean -a \
    && echo "source activate plague-phylogeography-qc" > ~/.bashrc

ENV PATH /opt/conda/envs/plague-phylogeography-qc/bin:$PATH
