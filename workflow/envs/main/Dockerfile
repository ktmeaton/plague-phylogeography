FROM continuumio/miniconda3

ARG GIT_COMMIT=unspecified

LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for ktmeaton/plague-phylogeography:main." \
      git_commit=$GIT_COMMIT

COPY main.yaml .
RUN conda install -c conda-forge mamba \
    && mamba env create -f main.yaml \
    && mamba clean -a \
    && echo "source activate plague-phylogeography" > ~/.bashrc

ENV PATH /opt/conda/envs/plague-phylogeography/bin:$PATH
