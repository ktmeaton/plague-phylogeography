FROM continuumio/miniconda3

LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for ktmeaton/plague-phylogeography: eager."

COPY eager.yaml .
RUN conda install -c conda-forge mamba \
    && conda env create -f eager.yaml \
    && conda clean -a \
    && echo "source activate plague-phylogeography-eager" > ~/.bashrc

ENV PATH /opt/conda/envs/plague-phylogeography-eager/bin:$PATH
