#!/bin/bash

# A failsafe script to move a LOCAL folder to a GCS bucket using `gcloud storage mv` or `gcloud storage cp`,
# retrying until all files are moved. Uses `tree` to ensure that no files remain.

set -euo pipefail

# Check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <source-dir> <gcs-destination> <command>"
    exit 1
fi

SOURCE_DIR="$1"
DEST="$2"
COMMAND="${3:-mv}" # command is default mv, but it can be overridden with cp

# Command must be mv or cp!
if [[ "$COMMAND" != "mv" && "$COMMAND" != "cp" ]]; then
    echo "Error: Command must be 'mv' or 'cp'."
    exit 1
fi

# Check that source is a local directory
if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory '$SOURCE_DIR' does not exist or is not a directory."
    exit 1
fi

# Check if `tree` is installed
if ! command -v tree &>/dev/null; then
    echo "Error: 'tree' command not found. Please install it (e.g., 'sudo apt install tree')."
    exit 1
fi

MAX_RETRIES=10
RETRY_DELAY=5
attempt=1

# Get absolute path
SOURCE_DIR_ABS=$(realpath "$SOURCE_DIR")

echo "Starting move from local '$SOURCE_DIR_ABS' to GCS '$DEST'"

# Initial count. Last line of `tree` output is expected to be in format `<N> directories, <file_count> files`
file_count=$(tree "$SOURCE_DIR_ABS" | tail -n 1 | awk '{print $3}')

while [[ "$file_count" != "0" && "$attempt" -le "$MAX_RETRIES" ]]; do
    echo "Attempt $attempt: $file_count files remain. Running gcloud storage $COMMAND..."
    gcloud storage "$COMMAND" -r "$SOURCE_DIR_ABS" "$DEST" || true # Ignore errors to allow retries

    sleep "$RETRY_DELAY"

    file_count=$(tree "$SOURCE_DIR_ABS" | tail -n 1 | awk '{print $3}')
    echo "Remaining files after attempt $attempt: $file_count"
    ((attempt++))
done

if [[ "$file_count" == "0" ]]; then
    echo "SUCCES! All files successfully moved."
else
    echo "ERROR! Maximum retries exceeded. $file_count files still remain in '$SOURCE_DIR_ABS'."
    exit 2
fi
