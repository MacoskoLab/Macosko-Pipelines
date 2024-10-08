FROM debian:bookworm
SHELL ["/bin/bash", "-c"]

# Install system libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends zip unzip gzip less wget curl    && \
    apt-get install -y --no-install-recommends sudo tmux htop dstat socat       && \
    apt-get install -y --no-install-recommends procps moreutils gnupg           && \
    apt-get install -y --no-install-recommends gdebi-core build-essential alien && \
    apt-get install -y --no-install-recommends libhdf5-dev hdf5-tools libpng-dev libtiff5-dev libjpeg-dev && \
    apt-get install -y --no-install-recommends apt-transport-https ca-certificates libssl-dev libxml2-dev zlib1g-dev && \
    apt-get install -y --no-install-recommends libfreetype6-dev libcairo2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev

# Install Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && apt-get update -y && apt-get install google-cloud-sdk -y

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

# Install Julia
RUN wget -P /opt https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.5-linux-x86_64.tar.gz && \
    tar -xzvf /opt/julia-1.10.5-linux-x86_64.tar.gz -C /opt                                             && \
    rm /opt/julia-1.10.5-linux-x86_64.tar.gz                                                            && \
    ln -s /opt/julia-1.10.5/bin/julia /usr/local/bin/julia

# Install Julia libraries
RUN julia -e 'using Pkg;                     \
              Pkg.add(["CSV",                \
                       "GZip",               \
                       "SIMD",               \
                       "HDF5",               \
                       "Plots",              \
                       "FASTX",              \
                       "Peaks",              \
                       "IJulia",             \
                       "ArgParse",           \
                       "CodecZlib",          \
                       "PDFmerger",          \
                       "IterTools",          \
                       "StatsBase",          \
                       "DataFrames",         \
                       "StatsPlots",         \
                       "StringViews",        \
                       "Combinatorics",      \
                       "Distributions",      \
                       "LinearAlgebra",      \
                       "LoopVectorization"])'

# Install R
RUN export R_VERSION=4.3.3 && \
    curl -O https://cdn.rstudio.com/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb && \
    echo y | gdebi r-${R_VERSION}_1_amd64.deb                                   && \
    rm r-${R_VERSION}_1_amd64.deb                                               && \
    sudo ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R                       && \
    sudo ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# Install R libraries
RUN R -e "install.packages(c('qs',          \
                             'glue',        \
                             'qpdf',        \
                             'dplyr',       \
                             'purrr',       \
                             'furrr',       \
                             'rlist',       \
                             'rdist',       \
                             'tidyr',       \
                             'hdf5r',       \
                             'dbscan',      \
                             'future',      \
                             'Seurat',      \
                             'cowplot',     \
                             'ggplot2',     \
                             'ggrastr',     \
                             'stringr',     \
                             'viridis',     \
                             'jsonlite',    \
                             'magrittr',    \
                             'gridExtra',   \
                             'SeuratObject'), repos='http://cran.rstudio.com/')"
RUN R -e "if (!require('BiocManager',quietly=T)) {install.packages('BiocManager',repos='http://cran.rstudio.com/')}; \
          BiocManager::install('rhdf5')"

# Install python
# Python 3.11.2 already installed
RUN sudo ln -s /usr/bin/python3 /usr/bin/python

# Install micromamba
RUN "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
# RUN echo "micromamba activate base" >> ~/.bashrc

# Install python packages
RUN /root/.local/bin/micromamba install -c conda-forge jupyterlab
