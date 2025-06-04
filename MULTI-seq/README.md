### Step 1: filterfastq

```
qsub -o /broad/macosko/jilong/nicolas/filterfastq/run.log \
     -l h_vmem=90g -notify -l h_rt=16:0:0 -j y -P macosko_lab -l os=RedHat7 \
     /broad/macosko/jilong/nicolas/filterfastq/run.sh /broad/macosko/jilong/nicolas/filterfastq \
     <in_R1.fastq> <in_R2.fastq> \
     <barcodes.tsv> \
     <out_R1.fastq> <out_R2.fastq> \
     1 true true
```


### Step 2: CITE-seq-Count

```pip install CITE-seq-Count==1.4.5```

```CITE-seq-Count -R1 TAGS_R1.fastq.gz -R2 TAGS_R2.fastq.gz -t TAG_LIST.csv -cbf X1 -cbl X2 -umif Y1 -umil Y2 -cells EXPECTED_CELLS -o OUTFOLDER```


### Helpful links
[MULTI-seq Sample Multiplexing for Single Cell Analysis and Sequencing](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/genomics/sequencing/multi-seq-sample-multiplexing-single-cell-analysis-sequencing)

[CITE-seq-Count](https://github.com/Hoohm/CITE-seq-Count)

| Binary            | line 124       | line 171   |
| ----------------- | -------------- | ---------- |
| filter_fastq      | substr(0, 16)  | substr(16) |
| filter_fastq_v2   | substr(0, 41)  | substr(41) |
| filter_fastq_v3   | substr(0, 32)  | substr(32) |
| filter_fastq_PONI | substr(46, 59) | substr(16) |
