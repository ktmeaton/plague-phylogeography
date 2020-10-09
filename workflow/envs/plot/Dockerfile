FROM continuumio/miniconda3
LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for plague-phylogeography: plot."
COPY plot.yaml /
RUN conda install -c conda-forge mamba; \
    mamba env create -f /plot.yaml; \
    conda clean -a; \
    echo "source activate plague-phylogeography-plot" > ~/.bashrc;
