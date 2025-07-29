FROM nvidia/cuda:12.1.1-cudnn8-devel-ubuntu22.04
#    nvidia/cuda:12.9.1-cudnn-devel-ubuntu22.04

# Install system packages
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y  \
    sudo zip unzip less file tree expect  \
    wget curl rsync socat                 \
    vim nano tmux screen htop ncdu dstat  \
    procps moreutils gnupg ssh git-all    \
    iproute2 net-tools lsof

# Install Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && apt-get update -y && apt-get install google-cloud-sdk -y

# Install micromamba
RUN curl -L micro.mamba.pm/install.sh | /bin/bash
# RUN /root/.local/bin/micromamba clean --locks

# Install python packages
RUN /root/.local/bin/micromamba install -y -n base -c conda-forge  \
                                  python=3.10.11 jupyterlab spyder \
                                  numpy pandas scipy scikit-learn  \
                                  matplotlib seaborn plotly pypdf  \
                                  pillow scikit-image tifffile     \
                                  opencv open3d && \
    /root/.local/bin/micromamba install -y -n base -c pytorch -c nvidia \
                    pytorch torchvision torchaudio pytorch-cuda=12.1 && \
    /root/.local/bin/micromamba run pip install cellpose==3.1.1.2 && \
    /root/.local/bin/micromamba clean --all --yes
# pip uninstall torch

# Set default theme to JupyterLab Dark
RUN DIR='/root/.jupyter/lab/user-settings/@jupyterlab/apputils-extension' && mkdir -p "$DIR" && echo '{\n  "theme": "JupyterLab Dark"\n}' > "$DIR/themes.jupyterlab-settings"


# https://hub.docker.com/r/nvidia/cuda
# https://catalog.ngc.nvidia.com/orgs/nvidia/containers/cuda
# https://pytorch.org/get-started/previous-versions/
# https://github.com/MouseLand/cellpose

# podman build --format docker -f Dockerfile -t gpu-image
# podman tag gpu-image us-central1-docker.pkg.dev/velina-208320/terra/gpu-image
# podman push us-central1-docker.pkg.dev/velina-208320/terra/gpu-image

# use UGER
# qrsh -verbose -q broad -now n -P macosko_lab -l h_rt=100:00:00 -l h_vmem=500g -l os=RedHat8 -l gpu=1 -hard -w e
# # -pe smp 8 -binding linear:8
# eval $GET_PORT
# podman run --rm -it --init --pull never \
#            --device nvidia.com/gpu=all  \
#            --security-opt=label=disable \
#            -v /broad/macosko:/broad/macosko:ro \
#            -v /broad/macosko/$USER:/broad/macosko/$USER:rw \
#            -p ${PORT_NUM:-8787}:8787 \
#            -w /broad/macosko/$USER \
#            gpu-image
# micromamba run jupyter lab --allow-root --ip='*' --port='8787' --NotebookApp.token='' --no-browser