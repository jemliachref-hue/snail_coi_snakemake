FROM condaforge/miniforge3:latest

WORKDIR /pipeline

COPY environment.yml .
RUN mamba env create -f environment.yml

COPY . .

SHELL ["conda", "run", "-n", "snail-coi", "/bin/bash", "-c"]

CMD ["conda", "run", "-n", "snail-coi", "snakemake", "--cores", "1"]
