### Step 0: Decompress the 10x barcode inclusion list

```gunzip -c barcodes.tsv.gz > barcodes.tsv```

### Step 1: CITE-seq-Count
```
CITE-seq-Count -R1 SAMPLE_R1_001.fastq.gz \
               -R2 SAMPLE_R2_001.fastq.gz \
               -t /broad/macosko/Nicolas/tag__List_OL478_OL765.csv \
               -cbf 1 -cbl 16 -umif 17 -umil 28 \
               -cells $(wc -l < barcodes.tsv) \
               -wl barcodes.tsv \
               -T $(nproc)
```

### Step 2: Run multiseq_pipeline.R
```
Arguments:
-i: Path to CITE-seq-Count outputs
-b: Path to Barcode Triplet File. Default: ./NL9_1_1009_combo.csv
-s: Samplesheet encoding mapping of MULTI-seq group id to condition. [Not required] See example: example_mulitseq_mapping.csv


Outputs:


```

### Helpful links
[MULTI-seq Sample Multiplexing for Single Cell Analysis and Sequencing](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/genomics/sequencing/multi-seq-sample-multiplexing-single-cell-analysis-sequencing)

[CITE-seq-Count](https://github.com/Hoohm/CITE-seq-Count)
