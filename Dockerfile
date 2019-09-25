FROM continuumio/anaconda:5.2.0

RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  libgsl-dev \
  libncurses-dev \
  libxt6 \
  libz-dev

RUN conda install -y python=2.7
RUN conda install -y biopython pandas psutil scons seaborn zlib
RUN conda install -y -c bioconda pysam
RUN conda install -y -c biocore mafft
RUN pip install colored-traceback dendropy==3.12.3
RUN pip install seqmagick==0.6.2

COPY . /templatedmutagenesis1
WORKDIR /templatedmutagenesis1/partis
RUN ./bin/build.sh
RUN conda install -y r-essentials \
    && unset R_LIBS_SITE \
    && R --vanilla --slave -e \
    'install.packages(c("argparse", "broom.mixed", "cowplot", "lme4", "gridExtra"), \
    repos="http://cran.rstudio.com/")'
WORKDIR /templatedmutagenesis1/PyMotifFinder
RUN pip install .
WORKDIR /templatedmutagenesis1
