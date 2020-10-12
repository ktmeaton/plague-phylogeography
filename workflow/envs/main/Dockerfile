FROM continuumio/miniconda3

ARG GIT_COMMIT=unspecified

LABEL authors="Katherine Eaton" \
      description="Docker image containing dependencies for ktmeaton/plague-phylogeography:main." \
      git_commit=$GIT_COMMIT

COPY main.yaml .
RUN conda install -c conda-forge mamba \
    && mamba env create -f main.yaml \
    && mamba clean -a \
    && echo "source activate plague-phylogeography" > ~/.bashrc;

ENV PATH /opt/conda/envs/plague-phylogeography/bin:$PATH

# Capture conda env details
RUN conda env export --name plague-phylogeography > plague-phylogeography.yaml

# Configure sra toolkit
RUN apt-get update \
    && apt-get --quiet install --yes uuid-runtime \
    && apt-get clean \
    && printf '/LIBS/GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg;

# Configure nf-core/eager
RUN nextflow run nf-core/eager -r dev -profile test_tsv,conda --fasta 'https://github.com/nf-core/test-datasets/raw/eager/reference/Mammoth/Mammoth_MT_Krause.fasta.gz';
