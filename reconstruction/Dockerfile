# https://docs.rapids.ai/install#selector
FROM nvcr.io/nvidia/rapidsai/base:24.04-cuda12.2-py3.11
SHELL ["/bin/bash", "-c"]
CMD ["/bin/bash"]
USER root

# Install system libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends zip unzip gzip less socat wget curl && \
    apt-get install -y --no-install-recommends sudo procps dstat moreutils gnupg   && \
    apt-get install -y --no-install-recommends python3-pip libhdf5-dev zlib1g-dev  && \
    apt-get install -y --no-install-recommends gdebi-core build-essential alien    && \
    apt-get install -y --no-install-recommends libpng-dev libtiff5-dev libjpeg-dev && \
    apt-get install -y --no-install-recommends apt-transport-https ca-certificates libssl-dev libxml2-dev && \
    apt-get install -y --no-install-recommends libfreetype6-dev libcairo2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev

# fix dstat bug (https://bugs.launchpad.net/ubuntu/+source/dstat/+bug/1905665)
RUN sed -i '547s/.*/            if isinstance(self.val[name], (tuple, list)):/g' /usr/bin/dstat
RUN sed -i '552s/.*/            elif isinstance(self.val[name], str):/g' /usr/bin/dstat

# Install Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && apt-get update -y && apt-get install google-cloud-sdk -y

# Install additional python packges
RUN pip install argparse seaborn plotly editdistance mappy umi_tools
RUN pip install umap-learn pydiffmap datafold MDAnalysis # dbmap
