FROM quay.io/broadinstitute/macosko-pipelines_base:current

COPY reconstruction/*.jl reconstruction/*.py /opt/Macosko-Pipelines/reconstruction/

RUN ln -s /opt/Macosko-Pipelines/reconstruction /usr/local/bin/reconstruction
