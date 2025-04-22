FROM debian:buster
SHELL ["/bin/bash", "-c"]

# Install system libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends zip unzip gzip less socat wget curl && \
    apt-get install -y --no-install-recommends sudo procps dstat moreutils gnupg   && \
    apt-get install -y --no-install-recommends gdebi-core build-essential alien    && \
    apt-get install -y --no-install-recommends apt-transport-https ca-certificates

# Install Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update -y && apt-get install google-cloud-sdk -y

# Add and install BCL Convert (RPM)
ADD bcl-convert-4.2.7-2.el7.x86_64.rpm /software/
RUN alien -i /software/bcl-convert-4.2.7-2.el7.x86_64.rpm && \
    rm /software/bcl-convert-4.2.7-2.el7.x86_64.rpm

RUN bcl-convert --version

# Add Cell Ranger
ADD cellranger-8.0.1 /software/cellranger-8.0.1
ENV PATH=/software/cellranger-8.0.1:$PATH

# Add Cell Ranger ATAC
ADD cellranger-atac-2.1.0 /software/cellranger-atac-2.1.0
ENV PATH=/software/cellranger-atac-2.1.0:$PATH

# Add Cell Ranger ARC
ADD cellranger-arc-2.0.2 /software/cellranger-arc-2.0.2
ENV PATH=/software/cellranger-arc-2.0.2:$PATH
