ARG BASE_IMAGE
FROM $BASE_IMAGE

ADD *.py /Macosko-Pipelines/bcl-convert/

# ../gcloud_build.sh -f Dockerfile -b BASE_IMAGE=us.gcr.io/mccarroll-scrna-seq/macosko-pipelines:2590bf6_1746626057 -t us.gcr.io/mccarroll-scrna-seq/macosko-pipelines-bcl-convert:$(git describe --match='' --always --abbrev=7 --dirty)_$(date +%s) --region=us-central1
