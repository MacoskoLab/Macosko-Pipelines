FROM debian:bookworm
ENV LANG C.UTF-8
ENV LANGUAGE C.UTF-8
ENV LC_ALL C.UTF-8

# Install system libraries
RUN apt-get update && apt-get install -y         \
    sudo wget curl rsync zip unzip less tree     \
    vim nano tmux screen htop dstat socat expect \
    procps moreutils gnupg iproute2 ssh git-all  \
    libudunits2-dev libgdal-dev                  \
    libncurses5-dev libncursesw5-dev             \
    zlib1g-dev liblzma-dev libbz2-dev            \
    gdebi-core build-essential cmake pkg-config alien \
    libhdf5-dev hdf5-tools libpng-dev libtiff5-dev libjpeg-dev \
    apt-transport-https ca-certificates libssl-dev libxml2-dev \
    libfreetype6-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libcairo2-dev

# Install Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && apt-get update -y && apt-get install google-cloud-sdk -y

# Install python
# (Python 3.11.2 already installed)
RUN sudo ln -s /usr/bin/python3 /usr/bin/python

# Install Samtools
RUN wget -P /opt https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xjvf /opt/samtools-1.21.tar.bz2 -C /opt && \
    rm /opt/samtools-1.21.tar.bz2 && \
    cd /opt/samtools-1.21 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install
# Install Bcftools
RUN wget -P /opt https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
    tar -xjvf /opt/bcftools-1.21.tar.bz2 -C /opt && \
    rm /opt/bcftools-1.21.tar.bz2 && \
    cd /opt/bcftools-1.21 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install
# Install HTSlib
RUN wget -P /opt https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar -xjvf /opt/htslib-1.21.tar.bz2 -C /opt && \
    rm /opt/htslib-1.21.tar.bz2 && \
    cd /opt/htslib-1.21 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install
# Install BEDTools
RUN wget -P /opt https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz && \
    tar -xzvf /opt/bedtools-2.29.1.tar.gz -C /opt && \
    rm /opt/bedtools-2.29.1.tar.gz && \
    cd /opt/bedtools2 && \
    make && \
    ln -s /opt/bedtools2/bin/bedtools /usr/local/bin/bedtools
# Install MMseqs2
RUN git clone https://github.com/soedinglab/MMseqs2.git /opt/MMseqs2 && \
    cd /opt/MMseqs2 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. .. && \
    make && \
    make install && \
    ln -s /opt/MMseqs2/build/bin/mmseqs /usr/local/bin/mmseqs

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

# Install Julia
RUN wget -P /opt https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.5-linux-x86_64.tar.gz && \
    tar -xzvf /opt/julia-1.10.5-linux-x86_64.tar.gz -C /opt && \
    rm /opt/julia-1.10.5-linux-x86_64.tar.gz && \
    ln -s /opt/julia-1.10.5/bin/julia /usr/local/bin/julia

# Install Julia libraries
RUN julia -e 'using Pkg;                \
              Pkg.add(["CSV",           \
                       "XAM",           \
                       "GZip",          \
                       "HDF5",          \
                       "Plots",         \
                       "FASTX",         \
                       "Peaks",         \
                       "IJulia",        \
                       "ArgParse",      \
                       "CodecZlib",     \
                       "PDFmerger",     \
                       "IterTools",     \
                       "StatsBase",     \
                       "DataFrames",    \
                       "StatsPlots",    \
                       "StringViews",   \
                       "Combinatorics", \
                       "Distributions", \
                       "LinearAlgebra"])'

# Install R
RUN export R_VERSION=4.3.3 && \
    curl -O https://cdn.rstudio.com/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb && \
    echo y | gdebi r-${R_VERSION}_1_amd64.deb && \
    rm r-${R_VERSION}_1_amd64.deb && \
    sudo ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R && \
    sudo ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# Install RStudio
RUN wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2024.09.0-375-amd64.deb && \
    echo y | gdebi rstudio-server-2024.09.0-375-amd64.deb && \
    rm rstudio-server-2024.09.0-375-amd64.deb && \
    rstudio-server stop

# Install R libraries
RUN R -e "install.packages(c('tidyverse', \
                             'dplyr', 'tidyr', 'purrr', 'magrittr', \
                             'future', 'furrr', 'parallelly', \
                             'data.table', 'rlist', 'stringr', 'stringi', 'glue', \
                             'ggplot2', 'ggrastr', 'cowplot', 'gridExtra', \
                             'scCustomize', 'viridis', \
                             'Seurat', 'SeuratObject', \
                             'rdist', 'sf', 'dbscan', \
                             'jsonlite', 'hdf5r', 'qpdf', 'qs', 'qs2', \
                             'devtools', 'remotes', 'R.utils', \
                             'shiny', 'IRkernel'), \
                             repos='http://cloud.r-project.org', \
                             Ncpus=$(nproc)L)"
# Install Bioconductor packages
RUN R -e "if (!require('BiocManager', quietly=T)) {install.packages('BiocManager', repos='http://cloud.r-project.org')}; \
          BiocManager::install(c('IRanges', 'GenomicRanges', 'GenomicFeatures', 'GenomicAlignments', \
                                 'ShortRead', 'Rsamtools', 'VariantAnnotation', 'rtracklayer', \
                                 'AnnotationDbi', 'BiocParallel', 'rhdf5' \
                                ), Ncpus=$(nproc)L)"

# Install micromamba
RUN curl -L micro.mamba.pm/install.sh | /bin/bash

# Install python packages
RUN /bin/bash -lc "micromamba install -c conda-forge jupyterlab \
                   numpy pandas scipy scikit-learn \
                   matplotlib seaborn plotly pypdf \
                   networkx rustworkx igraph graph-tool \
                   pynndescent umap-learn leidenalg"

# Install IRkernel
RUN /bin/bash -lc "micromamba run R -e 'IRkernel::installspec(user = FALSE)'"

ENTRYPOINT ["/bin/bash", "-lc"]
CMD ["/bin/bash", "-i"]

# podman build -f Dockerfile -t pipeline-image .
# podman save -o /broad/macosko/discopipeline/scripts/pipeline-image.tar pipeline-image:latest
