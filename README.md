# Templated Mutagenesis 1

[![Docker Repository on Quay](https://quay.io/repository/matsengrp/templatedmutagenesis1/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/templatedmutagenesis1)

This repository contains all of the code required to reproduce the figures and tables in the paper _Lack of evidence for a substantial rate of templated mutagenesis in B cell diversification_ by Julia Fukuyama, Branden J Olson, and Frederick A Matsen IV.


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


## Supplementary analysis

To run the analysis discussed in section "Consistent results using filtered donor gene sets" of the paper, you will need to install several python and R packages as laid out in the "Scripts" section of `data/reference_sets/README.md`.
Then, run the following commands which make use of the same `data` folder obtained from the first two steps in the Main analysis above.


    cd data/reference_sets/scripts
    sh make_tree.sh run_supplement
    cd ../../..
    mkdir supplementary_output
    docker run --mount src=${PWD}/data,target=/templatedmutagenesis1/data,type=bind,readonly=true \
        --mount src=${PWD}/supplementary_output,target=/templatedmutagenesis1/supplementary_output,type=bind \
        -t matsengrp/templatedmutagenesis1 scons MAKE_PLOTS=false OUTPUT_DIR=supplementary_output

This analysis takes several hours to complete.
Once finished, the results should be available in a `supplementary_output` directory.
Note that during the above, we made use of non-deterministic functions in `pyvolve` without the ability to provide the random seed, so your results may differ from those shown in Supplementary Tables III and IV.
Despite this caveat, the results should still look highly similar, and the conclusions should still hold.

To recover the original unfiltered dataset, you will need to rerun the `tar xzf templated-mutagenesis-data-v3.tar.gz` from above.
