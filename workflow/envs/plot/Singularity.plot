Bootstrap: docker

From: continuumio/miniconda3

%files
    plot.yaml

%post
    /opt/conda/bin/conda install -c conda-forge mamba
    mamba env create -f plot.yaml

%runscript
    exec /opt/conda/envs/$(head -n 1 plot.yaml | cut -f 2 -d ' ')/bin/"$@"
