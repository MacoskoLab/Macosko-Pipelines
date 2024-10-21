#!/bin/bash

# Read the command line
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <BCL> <INDEX> <TRANSCRIPTOME> <CHEMISTRY>"
    exit 1
fi

# Assert no spaces
for arg in "$@"; do
    if [[ "$arg" =~ " " ]]; then
        echo "Error: Argument '$arg' contains spaces."
        exit 1
    fi
done

BCL="$1"
INDEX="$2"
TRANSCRIPTOME="$3"
CHEMISTRY="$4"

ROOT="/broad/macosko/data/discopipeline"
BINARY="/broad/macosko/data/discopipeline/software/cellranger-8.0.1/bin/cellranger"

SBATCH_PARAMS="-C=RedHat7 -o=$ROOT/logs/$BCL/$INDEX/cellranger-count.out -J=cellranger-count-$INDEX -c=32 --mem=96G -t 96:00:00"
CELLRANGER_PARAMS="--id=$INDEX \
                   --output-dir=$ROOT/cellranger-count/$BCL \
                   --transcriptome=$ROOT/references/$TRANSCRIPTOME \
                   --fastqs=$ROOT/fastqs/$BCL \
                   --sample=$INDEX \
                   --chemistry=$CHEMISTRY \
                   --create-bam=true --include-introns=true --nosecondary"

sbatch $SBATCH_PARAMS $BINARY count $CELLRANGER_PARAMS
