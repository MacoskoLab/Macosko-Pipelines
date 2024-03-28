FROM debian:buster
SHELL ["/bin/bash", "-c"]

# Install system libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends dstat moreutils zip unzip gzip less socat wget curl && \
    apt-get install -y --no-install-recommends sudo gdebi-core build-essential alien && \
    apt-get install -y --no-install-recommends apt-transport-https ca-certificates gnupg

# Install Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && apt-get update -y && apt-get install google-cloud-sdk -y

# Add bcl2fastq2
ADD bcl2fastq2-v2-20-0-linux-x86-64.zip /software/
RUN unzip -d /software/ /software/bcl2fastq2-v2-20-0-linux-x86-64.zip && \
    alien -i /software/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && \
    rm /software/bcl2fastq2-v2*

# Add Cell Ranger
ADD cellranger-8.0.0 /software/cellranger-8.0.0
ENV PATH=/software/cellranger-8.0.0:$PATH
ADD cellranger-atac-2.1.0 /software/cellranger-atac-2.1.0
ENV PATH=/software/cellranger-arc-2.0.2:$PATH
ADD cellranger-arc-2.0.2 /software/cellranger-arc-2.0.2
ENV PATH=/software/cellranger-atac-2.1.0:$PATH
