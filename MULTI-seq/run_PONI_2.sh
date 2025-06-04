#!/bin/bash

source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3
reuse Java-1.8
source activate slideseq_pipeline_env

submission=$0
script_folder=$1
read1=$2
read2=$3
barcode=$4
outfile1=$5
outfile2=$6
distance_threshold=$7
correctbarcode=$8
uniquematch=$9

echo ${submission}

if [[ ${distance_threshold} == 0 ]]; then
	${script_folder}/filter_fastq_exactmatch ${read1} ${read2} ${barcode} ${outfile1} ${outfile2} ${distance_threshold} ${correctbarcode} ${uniquematch}
else
	${script_folder}/filter_fastq_PONI_2 ${read1} ${read2} ${barcode} ${outfile1} ${outfile2} ${distance_threshold} ${correctbarcode} ${uniquematch}
fi

