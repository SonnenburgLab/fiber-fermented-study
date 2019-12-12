# Base Image
FROM continuumio/miniconda3:4.6.14
USER root:root

# Data location
RUN mkdir -p /mnt
WORKDIR /mnt

# Install system-level programs
RUN apt-get update && apt-get install -qqy --no-install-suggests --no-install-recommends bc pigz parallel make gcc build-essential zlib1g-dev libbz2-dev

# Configuring conda channels
RUN conda config --add channels defaults && conda config --add channels bioconda && conda config --add channels conda-forge

# Installing conda packages
RUN conda install -y -c conda-forge awscli capnproto && conda install -y -c bioconda fastqc bbmap mash snakemake bowtie2 && conda clean -ay 

# Install latest version of sourmash
RUN git clone https://github.com/dib-lab/sourmash.git && cd sourmash && pip install -r requirements.txt && python setup.py install && cd ../


RUN chown -R root * && chgrp -R root * && chmod -R +rx ./
COPY . .
RUN chmod +rx *

# Metadata
LABEL container.maintainer="Hannah Wastyk <hwastyk@stanford.edu>" \
      container.base.image="continuumio/miniconda3:4.6.14" \
      container.version="0.0.1" \
      software.name="fastqc, bbmap, mash, snakemake, bowtie2, sourmash"
