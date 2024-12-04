#!/bin/bash

# Read the command line
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <BCL-FULL-PATH>"
    exit 1
fi

BCLPATH="$1"               ; echo "BCL-FULL-PATH: $BCLPATH"
BCL=$(basename "$BCLPATH") ; echo "BCL: $BCL"

# Assert no whitespace
for arg in "$@"; do
    if [[ "$arg" =~ [[:space:]] ]]; then
        echo "ERROR: Argument '$arg' contains whitespace"
        exit 1
    fi
done

ROOT="/broad/macosko/discopipeline"
BINARY="/broad/macosko/discopipeline/software/bcl-convert-4.3.6-2.el7.x86_64/bin/bcl-convert"
LOGDIR="$ROOT/logs/$BCL"

BCLCONVERT_PARAMS="--bcl-input-directory=$BCLPATH \
                   --output-directory=$ROOT/fastqs/$BCL \
                   --sample-sheet=$ROOT/samplesheets/$BCL/SampleSheet.csv \
                   --strict-mode=true"

SBATCH_PARAMS="-C RedHat7 -o $LOGDIR/bcl-convert.log -J bcl-convert-$BCL \
               -c 32 --mem 128G --time 72:00:00 \
               --mail-user mshabet@broadinstitute.org --mail-type END,FAIL,REQUEUE,INVALID_DEPEND,STAGE_OUT,TIME_LIMIT"

mkdir -p $LOGDIR
cd $ROOT/fastqs
sbatch $SBATCH_PARAMS --wrap "$BINARY $BCLCONVERT_PARAMS"
