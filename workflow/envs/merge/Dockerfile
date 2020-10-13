FROM continuumio/miniconda3

ARG GIT_COMMIT=unspecified

LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for ktmeaton/plague-phylogeography:latest." \
      git_commit=$GIT_COMMIT

# Add conda environments
COPY environment.yaml .

# Configure sra toolkit
RUN mkdir -p /root/.ncbi/ \
    && apt-get update \
    && apt-get --quiet install --yes uuid-runtime \
    && apt-get clean \
    && printf '/LIBS/GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg;

# Create conda environment
RUN conda install -c conda-forge mamba \
    && mamba env create -f environment.yaml \
    && mamba clean -a \
    && echo "source activate plague-phylogeography" > ~/.bashrc;

ENV PATH /opt/conda/envs/plague-phylogeography/bin:$PATH

# Capture conda env details
RUN conda env export --name plague-phylogeography > plague-phylogeography.yaml
