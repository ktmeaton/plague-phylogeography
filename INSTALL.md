# Install website

## Install hugo extended

```bash
wget https://github.com/gohugoio/hugo/releases/download/v0.75.1/hugo_extended_0.75.1_Linux-64bit.tar.gz
tar -xvf hugo_extended_0.75.1_Linux-64bit.tar.gz
```

## Install golang

```bash
conda env create -f environment.yaml
conda activate academic
```

## Run local server

```bash
./hugo server
```
