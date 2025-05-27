### Step 1: filterfastq

```
qsub -o /broad/macosko/jilong/nicolas/filterfastq/run.log -l h_vmem=90g -notify -l h_rt=16:0:0 -j y -P macosko_lab -l os=RedHat7 /broad/macosko/jilong/nicolas/filterfastq/run.sh /broad/macosko/jilong/nicolas/filterfastq <Read1_of_multiseq_barcode> <Read2_of_multiseq_barcode> <Cell_barcodes_from_cellranger> <Output_read1> <Output read2> 1 true true
```

### Step 2: CITE-seq-Count

https://github.com/Hoohm/CITE-seq-Count
