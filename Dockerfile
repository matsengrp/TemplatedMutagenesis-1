FROM continuumio/anaconda

RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  libgsl-dev \
  libncurses-dev \
  libxt6 \
  libz-dev

RUN conda update -y conda
RUN conda install -y python=2.7
RUN conda install -y biopython pandas psutil scons seaborn zlib
RUN conda install -y -c bioconda pysam
RUN conda install -y -c biocore mafft
RUN pip install colored-traceback dendropy==3.12.3
RUN pip install seqmagick==0.6.2

COPY . /gcgcgc
WORKDIR /gcgcgc/partis
RUN ./bin/build.sh
RUN conda install -y r-essentials \
    && unset R_LIBS_SITE \
    && R --vanilla --slave -e \
    'install.packages(c("argparse", "cowplot", "lme4", "gridExtra", "phytools"), \
    repos="http://cran.rstudio.com/")'
WORKDIR /gcgcgc/PyMotifFinder
RUN pip install .
WORKDIR /gcgcgc