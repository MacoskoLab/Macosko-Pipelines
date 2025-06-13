#!/bin/bash

source /broad/software/scripts/useuse
reuse UGER

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

ROOT="/broad/macosko/pipelines"
BINARY="/broad/macosko/pipelines/software/bcl-convert-4.3.6-2.el8.x86_64/bin/bcl-convert"
OUTDIR="$ROOT/fastqs/$BCL"
LOGDIR="$ROOT/logs/$BCL"

BCLCONVERT_PARAMS="--bcl-input-directory=$BCLPATH \
                   --output-directory=$OUTDIR \
                   --sample-sheet=$ROOT/samplesheets/$BCL/SampleSheet.csv \
                   --strict-mode=true"

UGER_PARAMS="-l h_rt=48:00:00 -l h_vmem=16G -pe smp 8 -binding linear:8 \
             -l os=RedHat8 -N bcl-convert-$BCL -o $LOGDIR/bcl-convert.log -j y \
             -m eas -M macosko-pipelines@broadinstitute.org -P macosko_lab -w e -notify"

mkdir -p $LOGDIR
cd $ROOT/fastqs
qsub $UGER_PARAMS -b y $BINARY $BCLCONVERT_PARAMS
