# Templated Mutagenesis 1

[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/matsengrp/templatedmutagenesis1.svg)](https://hub.docker.com/r/matsengrp/templatedmutagenesis1)

This repo contains all of the code required to reproduce all the figures and tables in the conversion paper.

## How to use it

A Docker image containing partis, PyMotifFinder, their dependencies, and the R packages required to run the analysis scripts is available on Docker Hub, and can be obtained by running

    docker pull matsengrp/templatedmutagenesis1

That docker image is built automatically from this repository.
If you wish can also build it yourself via

    docker build -t matsengrp/templatedmutagenesis1 .

The following commands will reproduce the figures in the paper, placing them in an `output` directory:

    wget https://zenodo.org/record/3572361/files/templated-mutagenesis-data-v2.tar.gz
    tar xzf templated-mutagenesis-data-v2.tar.gz
    mkdir output
    docker run --mount src=${PWD}/data,target=/templatedmutagenesis1/data,type=bind,readonly=true \
        --mount src=${PWD}/output,target=/templatedmutagenesis1/output,type=bind \
        -t matsengrp/templatedmutagenesis1 scons

The analysis takes several hours to complete.
