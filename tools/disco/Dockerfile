FROM debian:bookworm

# Install system libraries
RUN apt-get update && apt-get install -y \
    sudo wget curl rsync zip unzip less tree    \
    vim nano tmux screen htop dstat socat       \
    procps moreutils gnupg ssh git-all          \
    libudunits2-dev libgdal-dev                 \
    gdebi-core build-essential pkg-config alien \
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

# Install python
# Python 3.11.2 already installed
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install micromamba
RUN curl -L micro.mamba.pm/install.sh | /bin/bash

# Install python packages
RUN /bin/bash -c "source ~/.profile ; micromamba activate base ; \
                  micromamba install -c conda-forge jupyterlab \
                 "

# Ensure the container doesn't exit immediately
# CMD ["tail", "-f", "/dev/null"]
