FROM python:3.11-slim

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        curl \
        git \
        build-essential \
        libhdf5-dev \
        pkg-config \
        apt-transport-https \
        ca-certificates \
        gnupg \
    && rm -rf /var/lib/apt/lists/*

# Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] \
    https://packages.cloud.google.com/apt cloud-sdk main" \
    > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update && apt-get install -y google-cloud-cli \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies — pin versions for reproducibility
RUN pip install --no-cache-dir \
    "cnmf==1.6" \
    "anndata>=0.11" \
    "scanpy>=1.10" \
    "scipy>=1.11" \
    "pandas>=2.0" \
    "numpy>=1.24,<2.0" \
    "scikit-learn>=1.3" \
    "matplotlib>=3.7" \
    "seaborn>=0.12" \
    "h5py>=3.9"

COPY scripts/ /scripts/
RUN chmod +x /scripts/*.py

WORKDIR /workspace
