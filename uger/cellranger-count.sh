#!/bin/bash

# Read the command line
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <BCL> <INDEX> <TRANSCRIPTOME> <CHEMISTRY>"
    exit 1
fi

BCL="$1"           ; echo "BCL: $BCL"
INDEX="$2"         ; echo "Index: $INDEX"
TRANSCRIPTOME="$3" ; echo "Transcriptome: $TRANSCRIPTOME"
CHEMISTRY="$4"     ; echo "Chemistry: $CHEMISTRY"

# Assert no whitespace
for arg in "$@"; do
    if [[ "$arg" =~ [[:space:]] ]]; then
        echo "ERROR: Argument '$arg' contains whitespace"
        exit 1
    fi
done

ROOT="/broad/macosko/discopipeline"
BINARY="/broad/macosko/discopipeline/software/cellranger-8.0.1/bin/cellranger"
OUTDIR="$ROOT/cellranger-count/$BCL"
LOGDIR="$ROOT/logs/$BCL/$INDEX"

CELLRANGER_PARAMS="--id $INDEX \
                   --output-dir "$OUTDIR/$INDEX" \
                   --transcriptome $ROOT/references/$TRANSCRIPTOME \
                   --fastqs $ROOT/fastqs/$BCL \
                   --sample $INDEX \
                   --chemistry $CHEMISTRY \
                   --create-bam true \
                   --include-introns true \
                   --nosecondary \
                   --disable-ui"

SBATCH_PARAMS="-C RedHat7 -o $LOGDIR/cellranger-count.log -J cellranger-count-$BCL-$INDEX \
               -c 16 --mem 128G --time 72:00:00 \
               --mail-user mshabet@broadinstitute.org --mail-type END,FAIL,REQUEUE,INVALID_DEPEND,STAGE_OUT,TIME_LIMIT"

mkdir -p $OUTDIR
mkdir -p $LOGDIR
cd $OUTDIR
sbatch $SBATCH_PARAMS --wrap "$BINARY count $CELLRANGER_PARAMS"
