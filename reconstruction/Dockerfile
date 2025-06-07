ARG BASE_IMAGE
FROM $BASE_IMAGE

ADD *.jl *.py /Macosko-Pipelines/reconstruction/

# ../gcloud_build.sh -f Dockerfile -b BASE_IMAGE=us.gcr.io/mccarroll-scrna-seq/macosko-pipelines:2590bf6_1746626057 -t us.gcr.io/mccarroll-scrna-seq/macosko-pipelines-reconstruction:$(git describe --match='' --always --abbrev=7 --dirty)_$(date +%s) --region=us-central1
