FROM almalinux:8

# Install gcloud CLI
RUN cd /opt && \
    curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-linux-x86_64.tar.gz && \
    tar -xf google-cloud-cli-linux-x86_64.tar.gz && \
    rm google-cloud-cli-linux-x86_64.tar.gz && \
    ./google-cloud-sdk/install.sh --path-update true

# Install bcl-convert
COPY bcl-convert-4.3.13-2.el8.x86_64.rpm /opt/
RUN dnf install -y /opt/bcl-convert-4.3.13-2.el8.x86_64.rpm && \
    dnf clean all && \
    rm /opt/bcl-convert-4.3.13-2.el8.x86_64.rpm

# Install Cell Ranger
COPY cellranger-9.0.1.tar.gz /opt/
RUN tar -xzf /opt/cellranger-9.0.1.tar.gz -C /opt && \
    rm /opt/cellranger-9.0.1.tar.gz && \
    ln -s /opt/cellranger-9.0.1/bin/cellranger /usr/local/bin/cellranger

# Install Cell Ranger ARC
COPY cellranger-arc-2.0.2.tar.gz /opt/
RUN tar -xzf /opt/cellranger-arc-2.0.2.tar.gz -C /opt && \
    rm /opt/cellranger-arc-2.0.2.tar.gz && \
    ln -s /opt/cellranger-arc-2.0.2/bin/cellranger-arc /usr/local/bin/cellranger-arc
    
# Install Cell Ranger ATAC
COPY cellranger-atac-2.2.0.tar.gz /opt/
RUN tar -xzf /opt/cellranger-atac-2.2.0.tar.gz -C /opt && \
    rm /opt/cellranger-atac-2.2.0.tar.gz && \
    ln -s /opt/cellranger-atac-2.2.0/bin/cellranger-atac /usr/local/bin/cellranger-atac

# use UGER
# /broad/software/scripts/ish -l h_vmem=16G -l h_rt=24:00:00 -l os='RedHat8'
# cd /broad/macosko/pipelines/software
# wget -O Dockerfile https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/bcl-convert/Dockerfile
# podman build -f Dockerfile -t software-image
# podman tag software-image us-central1-docker.pkg.dev/velina-208320/terra/software-image
# podman push us-central1-docker.pkg.dev/velina-208320/terra/software-image
