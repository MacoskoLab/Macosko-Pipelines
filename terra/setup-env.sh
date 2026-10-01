#!/usr/bin/env bash
# Creates a micromamba environment with the libraries needed to run
# terra/submit.py and terra/results.py.
set -euo pipefail

ENV_NAME="${1:-terra-tools}"
PYTHON_VERSION="3.11"

micromamba create -n "$ENV_NAME" -y python="$PYTHON_VERSION" -c conda-forge

micromamba run -n "$ENV_NAME" pip install \
    gspread \
    gspread-dataframe \
    pandas \
    openpyxl \
    google-cloud-storage \
    google-auth \
    firecloud

echo ""
echo "Environment '$ENV_NAME' created. Activate it with:"
echo "  eval \"\$(micromamba shell hook --shell bash)\""
echo "  micromamba activate $ENV_NAME"


# One-time auth setup:
#   1. gcloud auth login
#   2. gcloud auth application-default login
#   3. gcloud config set project <GCP-PROJECT>
#   4. Check roles/iam.serviceAccountTokenCreator on the SA
#   5. (submit.py only) Register account at app.terra.bio and join the target workspace
