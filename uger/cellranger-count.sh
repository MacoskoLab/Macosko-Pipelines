#!/bin/bash

source /broad/software/scripts/useuse
reuse UGER

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

ROOT="/broad/macosko/pipelines"
BINARY="/broad/macosko/pipelines/software/cellranger-8.0.1/bin/cellranger"
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

UGER_PARAMS="-l h_rt=96:00:00 -l h_vmem=16G -pe smp 8 -binding linear:8 \
             -l os=RedHat8 -N cellranger-count-$BCL-$INDEX -o $LOGDIR/cellranger-count.log -j y \
             -m eas -M macosko-pipelines@broadinstitute.org -P macosko_lab -w e -notify"

mkdir -p $OUTDIR
mkdir -p $LOGDIR
cd $OUTDIR
qsub $UGER_PARAMS -b y $BINARY count $CELLRANGER_PARAMS
