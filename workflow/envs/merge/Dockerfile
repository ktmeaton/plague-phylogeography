FROM continuumio/miniconda3

ARG GIT_COMMIT=unspecified

LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for ktmeaton/plague-phylogeography:latest." \
      git_commit=$GIT_COMMIT

# Add conda environments
COPY environment.yaml .

# Configure sra toolkit
RUN apt-get update \
    && apt-get --quiet install --yes procps \
    && apt-get clean;

# Create conda environment
RUN conda install -c conda-forge mamba \
    && mamba env create -f environment.yaml \
    && mamba clean -a \
    && echo "source activate plague-phylogeography" > ~/.bashrc;

ENV PATH /opt/conda/envs/plague-phylogeography/bin:$PATH
ENV IN_DOCKER_CONTAINER Yes

# Install nf-core eager pipeline
RUN nextflow pull nf-core/eager -r dev;

# Install nodejs modules
RUN npm install --global auspice@2.22.0

# Capture conda env details
RUN conda env export --name plague-phylogeography > plague-phylogeography.yaml
