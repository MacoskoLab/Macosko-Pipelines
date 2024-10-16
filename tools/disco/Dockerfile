FROM debian:bookworm

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

# Install R
RUN export R_VERSION=4.3.3 && \
    curl -O https://cdn.rstudio.com/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb && \
    echo y | gdebi r-${R_VERSION}_1_amd64.deb                                   && \
    rm r-${R_VERSION}_1_amd64.deb                                               && \
    sudo ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R                       && \
    sudo ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# Install RStudio
RUN wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2024.09.0-375-amd64.deb && \
    sudo gdebi rstudio-server-2024.09.0-375-amd64.deb && \
    rm rstudio-server-2024.09.0-375-amd64.deb

# Install python
# Python 3.11.2 already installed
RUN sudo ln -s /usr/bin/python3 /usr/bin/python

# Install micromamba
RUN "${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# Install jupyterlab
RUN /root/.local/bin/micromamba install -c conda-forge jupyterlab

ENTRYPOINT ["/bin/bash"]
