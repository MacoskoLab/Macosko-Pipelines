#!/bin/bash

# Read the command line
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <BCL> <INDEX> <TRANSCRIPTOME> <CHEMISTRY>"
    exit 1
fi

BCL="$1"           ; echo $BCL
INDEX="$2"         ; echo $INDEX
TRANSCRIPTOME="$3" ; echo $TRANSCRIPTOME
CHEMISTRY="$4"     ; echo $CHEMISTRY

# Assert no whitespace
for arg in "$@"; do
    if [[ "$arg" =~ [[:space:]] ]]; then
        echo "ERROR: Argument '$arg' contains whitespace"
        exit 1
    fi
done

ROOT="/broad/macosko/data/discopipeline"
BINARY="/broad/macosko/data/discopipeline/software/cellranger-8.0.1/bin/cellranger"

SBATCH_PARAMS="-C RedHat7 -o $ROOT/logs/$BCL/$INDEX/cellranger-count.log -J cellranger-count-$BCL-$INDEX \
               -c 32 --mem 96G --time 96:00:00 \
               --mail-user mshabet@broadinstitute.org --mail-type END,FAIL,REQUEUE,INVALID_DEPEND,STAGE_OUT,TIME_LIMIT"
CELLRANGER_PARAMS="--id $INDEX \
                   --output-dir $ROOT/cellranger-count/$BCL \
                   --transcriptome $ROOT/references/$TRANSCRIPTOME \
                   --fastqs $ROOT/fastqs/$BCL \
                   --sample $INDEX \
                   --chemistry $CHEMISTRY \
                   --create-bam true \
                   --include-introns true \
                   --nosecondary"

sbatch $SBATCH_PARAMS $BINARY count $CELLRANGER_PARAMS
