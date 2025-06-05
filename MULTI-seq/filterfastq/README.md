This step has been deprecated in favor of the CITE-seq-Count option `-wl WHITELIST, --whitelist WHITELIST`

```
qsub -o /broad/macosko/jilong/nicolas/filterfastq/run.log \
     -l h_vmem=90g -notify -l h_rt=16:0:0 -j y -P macosko_lab -l os=RedHat7 \
     /broad/macosko/jilong/nicolas/filterfastq/run.sh /broad/macosko/jilong/nicolas/filterfastq \
     <in_R1.fastq> <in_R2.fastq> \
     <barcodes.tsv> \
     <out_R1.fastq> <out_R2.fastq> \
     1 true true
```

| Binary            | line 124       | line 171   |
| ----------------- | -------------- | ---------- |
| filter_fastq      | substr(0, 16)  | substr(16) |
| filter_fastq_v2   | substr(0, 41)  | substr(41) |
| filter_fastq_v3   | substr(0, 32)  | substr(32) |
| filter_fastq_PONI | substr(46, 59) | substr(16) |
