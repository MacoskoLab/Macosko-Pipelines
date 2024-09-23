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
RUN pip install numpy pandas matplotlib scipy scikit-learn seaborn
RUN pip install networkx igraph rustworkx
RUN pip install gcsfs plotly editdistance mappy umi_tools PyPDF2
RUN pip install umap-learn[tbb] pydiffmap datafold MDAnalysis # dbmap
RUN pip install pynndescent scikit-image datashader bokeh holoviews colorcet
RUN pip install tensorflow annoy tqdm

# Install Julia
RUN wget -P /software https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz && \
    tar -xzvf /software/julia-1.8.5-linux-x86_64.tar.gz -C /software                                       && \
    rm /software/julia-1.8.5-linux-x86_64.tar.gz
RUN /software/julia-1.8.5/bin/julia -e 'using Pkg;               \
                                        Pkg.add("CSV");          \
                                        Pkg.add("FASTX");        \
                                        Pkg.add("Plots");        \
                                        Pkg.add("Peaks");        \
                                        Pkg.add("ArgParse");     \
                                        Pkg.add("CodecZlib");    \
                                        Pkg.add("PDFmerger");    \
                                        Pkg.add("IterTools");    \
                                        Pkg.add("StatsBase");    \
                                        Pkg.add("DataFrames");   \
                                        Pkg.add("StatsPlots");   \
                                        Pkg.add("StringViews");  \
                                        Pkg.add("LinearAlgebra");\
                                        Pkg.add("Distributions");\
                                        Pkg.add("Combinatorics")'
