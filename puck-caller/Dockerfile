FROM nvidia/cuda:12.9.1-cudnn-devel-ubuntu22.04

# Install system packages
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y  \
    sudo zip unzip less tree expect       \
    wget curl rsync socat                 \
    vim nano tmux screen htop ncdu dstat  \
    procps moreutils gnupg ssh git-all    \
    iproute2 net-tools lsof

# Install Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && apt-get update -y && apt-get install google-cloud-sdk -y

# Install micromamba
RUN curl -L micro.mamba.pm/install.sh | /bin/bash
RUN /root/.local/bin/micromamba clean --locks

# Install python packages
RUN /root/.local/bin/micromamba install -y -n base -c conda-forge \
                                python=3.12.10 jupyterlab spyder  \
                                numpy pandas scipy scikit-learn   \
                                matplotlib seaborn plotly pypdf   \
                                pillow scikit-image tifffile      \
                                opencv open3d cellpose
# RUN /root/.local/bin/micromamba clean --all --yes

# https://hub.docker.com/r/nvidia/cuda
# https://catalog.ngc.nvidia.com/orgs/nvidia/containers/cuda
# podman build -f Dockerfile -t gpu-image # --format docker
# podman tag gpu-image us-central1-docker.pkg.dev/velina-208320/terra/gpu-image
# podman push us-central1-docker.pkg.dev/velina-208320/terra/gpu-image
