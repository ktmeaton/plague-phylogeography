FROM continuumio/miniconda3

LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for ktmeaton/plague-phylogeography: align."

COPY align.yaml .
RUN conda install -c conda-forge mamba \
    && conda env create -f align.yaml \
    && conda clean -a \
    && echo "source activate plague-phylogeography-align" > ~/.bashrc

ENV PATH /opt/conda/envs/plague-phylogeography-align/bin:$PATH
