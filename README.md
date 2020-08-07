# Templated Mutagenesis 1

[![Docker Repository on Quay](https://quay.io/repository/matsengrp/templatedmutagenesis1/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/templatedmutagenesis1)

This repository contains all of the code required to reproduce the figures and tables in the paper [_Lack of evidence for a substantial rate of templated mutagenesis in B cell diversification_](http://dx.doi.org/10.4049/jimmunol.2000092) by Julia Fukuyama, Branden J Olson, and Frederick A Matsen IV.


## How to use it

A Docker image containing partis, PyMotifFinder, their dependencies, and the R packages required to run the analysis scripts is available on Quay, and can be obtained by running the following commands:

    docker pull quay.io/matsengrp/templatedmutagenesis1
    docker tag quay.io/matsengrp/templatedmutagenesis1 matsengrp/templatedmutagenesis1

That docker image is built automatically from this repository.
If you wish can also build it yourself via

    docker build -t matsengrp/templatedmutagenesis1 .


## Main analysis

The following commands will reproduce the figures in the paper, placing them in an `output` directory:

    wget https://zenodo.org/record/3741809/files/templated-mutagenesis-data-v3.tar.gz
    tar xzf templated-mutagenesis-data-v3.tar.gz
    mkdir output
    docker run --mount src=${PWD}/data,target=/templatedmutagenesis1/data,type=bind,readonly=true \
        --mount src=${PWD}/output,target=/templatedmutagenesis1/output,type=bind \
        -t matsengrp/templatedmutagenesis1 scons

The analysis takes several hours to complete.

Instructions for running the supplementary analysis as well as results can be found [here](docs/supplement.md).
