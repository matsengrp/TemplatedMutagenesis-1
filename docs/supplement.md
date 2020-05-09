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
Note that during the above, we made use of non-deterministic functions in `pyvolve` without the ability to provide the random seed, so your results may differ from those shown in the tables below.
Despite this caveat, the results should still look highly similar, and the conclusions should still hold.

To recover the original unfiltered dataset, you will need to rerun the `tar xzf templated-mutagenesis-data-v3.tar.gz` from above.

## Results

Here we include example tables from running the above analysis, which are exactly those discussed in section "Consistent results using filtered donor gene sets" of the JI paper.
In all of the following tables, _k_ denotes tract length, PyPolyMF rate is the naive PyPolyMF estimate of the rate of templated mutagenesis, PyPolyMF FPR is the PyPolyMF false positive rate, UB denotes upper bound, and the number in parentheses denotes the assumed sensitivity (true positive rate) of PyPolyMF.

### Mice

Bounds for mice, with open reading frame (ORF) and pseudogene (P) sequence reads _included in_ the donor set:

| _k_ | PyPMF rate | PyPMF FPR | UB (1) | UB (.99) | UB (.95) | UB (.9)|
|---|------------|-----------|--------|----------|----------|--------|
|8|0.79|0.83|0|0|0|0|
|9|0.54|0.5|0.08|0.08|0.09|0.1|
|10|0.28|0.25|0.05|0.05|0.05|0.05|
|11|0.15|0.14|0.01|0.01|0.01|0.02|
|12|0.11|0.09|0.01|0.01|0.01|0.02|
|13|0.07|0.07|0|0|0|0|
|14|0.06|0.05|0.01|0.01|0.01|0.01|
 
Bounds for mice, with open reading frame (ORF) and pseudogene (P) sequence reads _excluded from_ the donor set:

| _k_ | PyPMF rate | PyPMF FPR | UB (1) | UB (.99) | UB (.95) | UB (.9)|
|---|------------|-----------|--------|----------|----------|--------|
|8|0.62|0.7|0|0|0|0|
|9|0.33|0.35|0|0|0|0|
|10|0.19|0.17|0.02|0.02|0.02|0.02|
|11|0.11|0.09|0.02|0.02|0.02|0.02|
|12|0.07|0.06|0.01|0.01|0.01|0.01|
|13|0.05|0.05|0|0|0|0|
|14|0.04|0.03|0.01|0.01|0.01|0.01|

### Humans

Bounds for humans, with open reading frame (ORF) and pseudogene (P) sequence reads _included in_ the donor set:

| _k_ | PyPMF rate | PyPMF FPR | UB (1) | UB (.99) | UB (.95) | UB (.9)|
|---|------------|-----------|--------|----------|----------|--------|
|8|0.73|0.78|0|0|0|0|
|9|0.44|0.43|0.02|0.02|0.02|0.02|
|10|0.26|0.2|0.07|0.08|0.08|0.09|
|11|0.18|0.11|0.09|0.09|0.09|0.1|
|12|0.15|0.07|0.09|0.1|0.1|0.11|
|13|0.15|0.05|0.1|0.1|0.11|0.11|
|14|0.13|0.03|0.1|0.11|0.11|0.12|
 
 
Bounds for humans, with open reading frame (ORF) and pseudogene (P) sequence reads _excluded from_ the donor set:

| _k_ | PyPMF rate | PyPMF FPR | UB (1) | UB (.99) | UB (.95) | UB (.9)|
|---|------------|-----------|--------|----------|----------|--------|
|8|0.42|0.46|0|0|0|0|
|9|0.25|0.2|0.07|0.07|0.07|0.08|
|10|0.17|0.09|0.09|0.09|0.09|0.1|
|11|0.14|0.05|0.09|0.09|0.1|0.1|
|12|0.13|0.04|0.09|0.09|0.1|0.1|
|13|0.13|0.03|0.1|0.1|0.1|0.11|
|14|0.12|0.02|0.1|0.1|0.11|0.11|
