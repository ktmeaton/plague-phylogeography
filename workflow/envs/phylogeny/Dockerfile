FROM continuumio/miniconda3

LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for ktmeaton/plague-phylogeography: phylogeny."

COPY phylogeny.yaml .
RUN conda install -c conda-forge mamba \
    && conda env create -f phylogeny.yaml \
    && conda clean -a \
    && echo "source activate plague-phylogeography-phylogeny" > ~/.bashrc

ENV PATH /opt/conda/envs/plague-phylogeography-phylogeny/bin:$PATH
