FROM debian:bookworm
ENV LANG C.UTF-8
ENV LANGUAGE C.UTF-8
ENV LC_ALL C.UTF-8

# Install system libraries
RUN apt-get update && apt-get install -y \
    sudo wget curl rsync zip unzip less tree     \
    vim nano tmux screen htop dstat socat expect \
    procps moreutils gnupg iproute2 ssh git-all  \
    libudunits2-dev libgdal-dev                  \
    gdebi-core build-essential pkg-config alien  \
    libhdf5-dev hdf5-tools libpng-dev libtiff5-dev libjpeg-dev \
    apt-transport-https ca-certificates libssl-dev libxml2-dev zlib1g-dev \
    libfreetype6-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libcairo2-dev

# Install Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && apt-get update -y && apt-get install google-cloud-sdk -y

# Install R
RUN export R_VERSION=4.3.3 && \
    curl -O https://cdn.rstudio.com/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb && \
    echo y | gdebi r-${R_VERSION}_1_amd64.deb                                   && \
    rm r-${R_VERSION}_1_amd64.deb                                               && \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R                            && \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# Install RStudio
RUN wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2024.09.0-375-amd64.deb && \
    echo y | gdebi rstudio-server-2024.09.0-375-amd64.deb && \
    rm rstudio-server-2024.09.0-375-amd64.deb && \
    rstudio-server stop

# Install R libraries
# RUN R -e "install.packages(c('qs',          \
#                              'glue',        \
#                              'qpdf',        \
#                              'dplyr',       \
#                              'purrr',       \
#                              'furrr',       \
#                              'rlist',       \
#                              'rdist',       \
#                              'tidyr',       \
#                              'hdf5r',       \
#                              'dbscan',      \
#                              'future',      \
#                              'Seurat',      \
#                              'cowplot',     \
#                              'ggplot2',     \
#                              'ggrastr',     \
#                              'stringr',     \
#                              'viridis',     \
#                              'jsonlite',    \
#                              'magrittr',    \
#                              'gridExtra',   \
#                              'SeuratObject'), repos='http://cran.rstudio.com/')"
# RUN R -e "if (!require('BiocManager',quietly=T)) {install.packages('BiocManager',repos='http://cran.rstudio.com/')}; \
#           BiocManager::install('rhdf5')"

# Install python
# Python 3.11.2 is already installed
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install micromamba
RUN curl -L micro.mamba.pm/install.sh | /bin/bash

# Install python packages
RUN /bin/bash -lc "micromamba install -c conda-forge jupyterlab \
                   numpy pandas scipy scikit-learn \
                   matplotlib seaborn plotly pypdf"

# Create password for RStudio and JupyterLab
ARG USER=root
ARG PASSWORD=password
# RStudio uses the linux password
RUN echo "$USER:$PASSWORD" | chpasswd
RUN echo "auth-minimum-user-id=0" >> /etc/rstudio/rserver.conf
# JupyterLab uses a config file
RUN expect -c 'spawn bash -lc "micromamba run jupyter lab password"; expect "Enter password:"; send "$env(PASSWORD)\r"; expect "Verify password:"; send "$env(PASSWORD)\r"; expect eof'

# Create bash functions
RUN printf "\nrstudio() { \n\
    tmux new -s rstudio -d bash -ic \"rstudio-server start ; exec bash\" \n\
    exec bash \n}\n" >> /root/.bashrc
RUN printf "\njupyterlab() { \n\
    tmux new -s jupyterlab -d exec bash -ic \"micromamba run jupyter lab --allow-root --ip='*' --port='8787' --NotebookApp.token='' --no-browser\" \n\
    exec bash \n}\n" >> /root/.bashrc

ENTRYPOINT ["/bin/bash", "-lc"]
CMD ["/bin/bash", "-i"]
