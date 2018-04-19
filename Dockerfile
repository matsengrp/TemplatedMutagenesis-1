FROM matsengrp/cpp

# ----------------------------------------------------------------------------------------
RUN sed -i 's/httpredir/ftp.us/' /etc/apt/sources.list
RUN apt-get update && apt-get install -y \
    libgsl0-dev \
    libncurses5-dev \
    libxml2-dev \
    libxslt1-dev \
    mafft \
    r-base
RUN pip install numpy  # putting them on different lines allows docker's caching to defeat pip's slowness
RUN pip install scipy
RUN pip install scikit-learn
RUN pip install matplotlib
RUN pip install pandas
RUN pip install biopython
RUN pip install dendropy==3.12.3
RUN pip install pysam
RUN pip install pyyaml
RUN pip install seaborn
RUN pip install colored_traceback
RUN pip install psutil
RUN pip install seqmagick==0.6.2

RUN R --vanilla --slave -e 'install.packages(c("TreeSim", "TreeSimGM", "bios2mds", "argparse", "reshape2", "ggplot2", "gridExtra", "phytools", "dplyr", "magrittr", "lme4", "broom"), repos="http://cran.rstudio.com/")'
# packages that are needed for packages we need to install from source
RUN R --vanilla --slave -e 'install.packages(c("data.table", "maps", "animation", "clusterGeneration", "msm", "numDeriv", "plotrix", "scatterplot3d", "quadprog", "igraph", "nnls", "fastmatch"), repos = "http://cran.rstudio.com")'
RUN R --vanilla --slave -e 'install.packages("http://cran.r-project.org/src/contrib/Archive/ape/ape_3.5.tar.gz", repos=NULL, type="source")'
RUN R --vanilla --slave -e 'install.packages("http://cran.r-project.org/src/contrib/Archive/phangorn/phangorn_2.0.4.tar.gz", repos=NULL, type="source")'
RUN R --vanilla --slave -e 'install.packages(c("http://cran.r-project.org/src/contrib/Archive/phytools/phytools_0.5-38.tar.gz"), repos=NULL, type = "source")'

COPY . /gcgcgc
WORKDIR /gcgcgc/partis
# build partis
RUN ./bin/build.sh
# install PyMotifFinder
WORKDIR /gcgcgc/PyMotifFinder
RUN pip install .
WORKDIR /gcgcgc
