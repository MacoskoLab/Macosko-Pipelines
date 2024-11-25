FROM debian:bookworm
ENV LANG C.UTF-8
ENV LANGUAGE C.UTF-8
ENV LC_ALL C.UTF-8
RUN echo "\numask 002" >> /etc/profile

# Install system libraries
RUN apt-get update && apt-get install -y         \
    sudo wget curl rsync zip unzip less tree     \
    vim nano tmux screen htop dstat socat expect \
    procps moreutils gnupg ssh git-all           \
    iproute2 net-tools lsof                      \
    libudunits2-dev libgdal-dev                  \
    libncurses5-dev libncursesw5-dev             \
    zlib1g-dev liblzma-dev libbz2-dev            \
    gdebi-core build-essential cmake pkg-config alien \
    libhdf5-dev hdf5-tools libpng-dev libtiff5-dev libjpeg-dev \
    apt-transport-https ca-certificates libssl-dev libxml2-dev \
    libfreetype6-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libcairo2-dev

# Install Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && apt-get update -y && apt-get install google-cloud-sdk -y

# Copy the key into the image
COPY velina-208320-bead39436991.json /root/velina-208320-key.json

# Authenticate
RUN gcloud auth activate-service-account --key-file=/root/velina-208320-key.json && \
    gcloud config set project velina-208320

# copied from fc-71ac3b81-2441-4171-8038-baf653634620
RUN mkdir -p ~/sc-online && \
    gcloud storage cp -r gs://macosko_data/sc-online/code ~/sc-online && \
    gcloud storage cp -r gs://macosko_data/sc-online/serverFiles ~/sc-online

# Install python
# (Python 3.11.2 already installed)
RUN sudo ln -s /usr/bin/python3 /usr/bin/python

# Install R (version in sc-online notebooks is 4.2.1 https://app.terra.bio/#workspaces/fbrihuman/sconline_integrative_analysis/analysis/launch/Settingup_sconline_source.ipynb)
RUN export R_VERSION=4.2.1 && \
    curl -O https://cdn.rstudio.com/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb && \
    echo y | gdebi r-${R_VERSION}_1_amd64.deb && \
    rm r-${R_VERSION}_1_amd64.deb && \
    sudo ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R && \
    sudo ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# Install RStudio
RUN wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2024.09.0-375-amd64.deb && \
    echo y | gdebi rstudio-server-2024.09.0-375-amd64.deb && \
    rm rstudio-server-2024.09.0-375-amd64.deb && \
    rstudio-server stop && sleep 1

# Install R libraries based on order of operations an packages found in 
# https://app.terra.bio/#workspaces/fbrihuman/sconline_integrative_analysis/analysis/launch/Settingup_sconline_source.ipynb
RUN R -e "install.packages('devtools', repos='http://cloud.r-project.org', Ncpus=$(nproc)L)"
RUN R -e "devtools::install_version('XML', version = '3.98-1.17', repos = 'http://cran.us.r-project.org', Ncpus=$(nproc)L)"


RUN R -e "install.packages(c( \
            'BiocManager', 'ggbeeswarm', 'CrossClustering', 'parallel', 'colorRamps', \
            'ggplot2', 'reshape', 'reshape2', 'Hmisc', 'patchwork', 'dplyr', \
            'ClusterR', 'hues', 'effsize', 'ggwordcloud'), \
        repos='http://cloud.r-project.org', \
        Ncpus=$(nproc)L)"

RUN R -e "BiocManager::install(c( \
            'Biobase', 'qvalue', 'lfa', 'multtest', 'limma', 'biomaRt', 'GenomicFeatures', \
            'EnsDb.Hsapiens.v75', 'monocle', 'org.Hs.eg.db', 'DropletUtils', 'scran', 'zoo', \
            'mixtools', 'rbokeh', 'DT', 'NMF', 'pheatmap', 'R2HTML', 'Rtsne', 'doMC', 'doRNG', \
            'EnsDb.Mmusculus.v79', 'SC3', 'PCAtools'), \
        Ncpus=$(nproc)L)"

RUN R -e "install.packages(c( \
            'Seurat', 'jackstraw', 'caret', 'IRkernel'), \
        repos='http://cloud.r-project.org', \
        Ncpus=$(nproc)L)"

RUN R -e "devtools::install_github('immunogenomics/harmony', Ncpus=$(nproc)L)"

RUN R -e "install.packages(c( \
            'ggraph', 'ape', 'seriation', 'igraph', 'qs', 'irlba', 'network', 'seriation', 'ggnetwork', \
            'cowplot', 'WGCNA', 'scatterpie'), \
        repos='http://cloud.r-project.org', \
        Ncpus=$(nproc)L)"

RUN R -e "BiocManager::install(c( \
            'impute', 'preprocessCore', 'GO.db', 'scater', 'EnsDb.Hsapiens.v86'), \
        Ncpus=$(nproc)L)"

# Install micromamba
RUN curl -L micro.mamba.pm/install.sh | /bin/bash
RUN /bin/bash -lc "micromamba clean --locks"

# Install python packages
RUN /bin/bash -lc "micromamba install -c conda-forge \
                    anndata jupyterlab leidenalg matplotlib \
                    numpy pandas plotly pypdf \
                    scikit-learn scipy seaborn umap-learn"

# Install IRkernel
RUN /bin/bash -lc "micromamba run R -e 'IRkernel::installspec(user = FALSE)'"

ENTRYPOINT ["/bin/bash", "-lc"]
CMD ["/bin/bash", "-i"]

# podman build -f Dockerfile -t pipeline-image .
# podman save -o /broad/macosko/discopipeline/scripts/pipeline-image.tar pipeline-image:latest