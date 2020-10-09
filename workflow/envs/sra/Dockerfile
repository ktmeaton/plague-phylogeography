FROM continuumio/miniconda3

LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for ktmeaton/plague-phylogeography: sra."

COPY sra.yaml .
RUN conda install -c conda-forge mamba \
    && conda env create -f sra.yaml \
    && conda clean -a \
    && echo "source activate plague-phylogeography-sra" > ~/.bashrc

ENV PATH /opt/conda/envs/plague-phylogeography-sra/bin:$PATH
