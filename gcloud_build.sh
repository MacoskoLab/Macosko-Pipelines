#!/bin/bash

# Enables passing build arguments to the Docker build command when using `gcloud builds submit`.
# Extended from https://til.simonwillison.net/cloudrun/using-build-args-with-cloud-run

# Example usage:
# ./gcloud_build.sh -f Dockerfile -b BASE_IMAGE=us.gcr.io/mccarroll-scrna-seq/macosko-pipelines:2590bf6_1746626057 -t us.gcr.io/mccarroll-scrna-seq/macosko-pipelines-reconstruction:$(git describe --match='' --always --abbrev=7 --dirty)_$(date +%s) --region=us-central1

set -euo pipefail

PROGNAME=$(basename "$0")

build_args=''
image_tag=''

function usage () {
    cat <<EOF
Usage: $PROGNAME [-h] -f Dockerfile -t image-tag  [-b build-arg]... [--] [gcloud-build-args]
Build and push a Docker image using Google Cloud.

Options:
  -h, --help:                   Show this help message
  -f, --file Dockerfile:        Path to the Dockerfile
  -t, --tag image-tag:          Tag for the built image
  -b, --build-arg NAME=VALUE:   Build argument for the Docker image
  gcloud-build-args:            Additional arguments to pass to \`gcloud builds submit\`
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        -f|--file)
            if [[ -z "$2" ]]; then
                echo "ERROR: Missing argument for $1" >&2
                usage
                exit 1
            fi
            docker_file="$2"
            build_args="$build_args, '--file', '$(basename "$2")'"
            shift 2
            ;;
        -t|--tag)
            if [[ -z "$2" ]]; then
                echo "ERROR: Missing argument for $1" >&2
                usage
                exit 1
            fi
            image_tag="$2"
            build_args="$build_args, '--tag', '$2'"
            shift 2
            ;;
        -b|--build-arg)
            if [[ -z "$2" ]]; then
                echo "ERROR: Missing argument for $1" >&2
                usage
                exit 1
            fi
            build_args="$build_args, '--build-arg', '$2'"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            break
    esac
done

if [[ -z "${docker_file:-}" ]]; then
    echo "ERROR: Dockerfile not specified" >&2
    usage
    exit 1
fi

if [[ -z "${image_tag:-}" ]]; then
    echo "ERROR: Image tag not specified" >&2
    usage
    exit 1
fi

cloudbuild_tmp_dir=$(mktemp -d /tmp/gcloud_build.XXXX)
cloudbuild_yaml="$cloudbuild_tmp_dir/cloudbuild.yaml"

function cleanup() {
    rm -rf "$cloudbuild_tmp_dir"
}

trap cleanup EXIT

builddir="$(dirname "$docker_file")"
rsync --archive --quiet "$builddir"/ "$cloudbuild_tmp_dir"/

cat <<EOF > "$cloudbuild_yaml"
steps:
- name: 'gcr.io/cloud-builders/docker'
  args: ['build'$build_args, '.']
- name: 'gcr.io/cloud-builders/docker'
  args: ['push', '$image_tag']
EOF

gcloud builds submit --config "$cloudbuild_yaml" "$@"
