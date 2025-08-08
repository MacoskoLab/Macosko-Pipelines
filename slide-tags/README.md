Quick-start commands
--------------------
```
julia spatial-count.jl fastq_dir puck_dir output
Rscript run-positioning.R gex_dir output/SBcounts.h5 output
```

Example input directory structure:
<pre>
slide-tags/
├── spatial-count.jl
├── run-positioning.R
├── positioning.R
├── helpers.R
├── plots.R
├── fastq_dir/
│   ├── SI-TT-A1_S1_L001_R1_001.fastq.gz
│   ├── SI-TT-A1_S1_L001_R2_001.fastq.gz
|   └── ...
├── puck_dir/
|   └── Puck.csv
└── gex_dir/
    ├── filtered_feature_bc_matrix.h5
    ├── molecule_info.h5
    └── metrics_summary.csv
</pre>
Output folder:
<pre>
output/
├── SBcounts.h5
├── matrix.csv.gz
├── spatial_metadata.json
├── coords.csv
├── seurat.qs
└── summary.pdf
</pre>

Seurat Metadata
---------------
* orig.ident: `Slide-tags`
* nCount_RNA: total number of molecules detected within a cell
* nFeature_RNA: number of genes detected in each cell
* logumi: `log10(nCount_RNA+1)`
* cb: cell barcode, with trailing `-[0-9]*` removed
* cb_index: 1-indexed cell barcode enumeration
* percent.mt: `PercentageFeatureSet` with `pattern="^(MT-|mt-)"`
* pct.intronic:
    * 10X: `(umi_type==0) / (umi_type==0 + umi_type==1)`
    * Optimus: `reads_mapped_intronic / (reads_mapped_intronic + reads_mapped_exonic)`
* RNA_snn_res.0.8: `FindNeighbors` with `dims=1:30`, then `FindClusters` with `resolution=0.8`
* seurat_clusters: same as above

Seurat Misc
-----------

SBcounts.h5
-----------
* lists
    * cb_list - cell barcode sequences
    * sb_list - spatial barcode sequences
    * puck_list - puck filenames
    * R1_list - R1 FASTQ filenames
    * R2_list - R2 FASTQ filenames
* matrix
    * cb_index - index of the cell barcode in lists/cb_list
    * umi - 2-bit encoded UMI sequence
    * sb_index - index of the cell barcode in lists/sb_list
    * reads - number of read sequences with this 3-tuple of barcodes
* puck
    * sb - spatial barcode sequence
    * x - x-coordinate of bead centroid
    * y - y-coordinate of bead centroid
    * puck_index - index into puck_list
* metadata
    * switch - False if R1 contains CB/UMI and R2 contains SB, True otherwise
    * bead - ID of read structure
    * num_lowQbeads - number of low-quality beads filtered from the puck
    * reads - total number of reads
    * reads_filtered - total number of reads passing the UMI, UP, and SB filters
    * R1/2_tooshort - number of reads with fewer than the allowable cycles
    * UMI - filtering statistics on the UMI sequence
    * UP - filtering statistics on the UP
    * SB - filtering statistics on the spatial barcode sequence
    * SB_HD - location of spatial barcode fuzzy match

Spatial Library Metadata
------------------------

* SB_info
    * R1s
    * R2s
    * pucks
    * UMI_downsampling
    * switch_R1R2
    * bead_type
    * remap_10X_CB
    * UMI_pct_in_called_cells
    * sequencing_saturation
* puck_info
    * puck_name - filename of puck
    * num_beads - total number of beads
    * num_lowQ - number of beads removed before FASTQ parsing
    * num_dup - number of beads removed by collision (puck multiplexing)
    * num_N - number of beads removed for having an N
    * num_degen - number of beads removed for being degnerate
    * puck_boundaries - x-coordinates demarcating puck start/end positions
    * umi_final - number of spatial UMIs associated with the puck after filtering
* UMI_filtering
    * N: number of reads where the UMI contained a N
    * homopolymer: number of reads where the UMI was degenerate
* UP_matching: exact, fuzzy, GG, none
* SB_matching: exact, HD1, HD1ambig, none
* SB_fuzzy_position: the number of barcodes which fuzzy-matched to each base
* CB_matching: exact, fuzzy, ambig, none
* SB_filtering
    * reads_total - total number of reads in the input FASTQ
    * reads_tooshort - FASTQ sequence not long enough
    * reads_noumi - UMI was low-quality
    * reads_noup - UP site not matched
    * reads_nosb - SB did not match puck whitelist
    * reads_lqsb - SB was low-quality
    * reads_nocb - CB did not match called cell whitelist
    * reads_final - number of reads after all above filters
    * UMIs_final - number of UMIs after all above filters


Bead structures
---------------
* `V10: JJJJJJJJ  TCTTCAGCGTTCCCGAGA JJJJJJJ  NNNNNNNVV`
* `V17: JJJJJJJJJ TCTTCAGCGTTCCCGAGA JJJJJJJJ NNNNNNNNN`
* `V15: JJJJJJJJJJJJJJJ   CTGTTTCCTG NNNNNNNNN`
* `V16: JJJJJJJJJJJJJJJJJ CTGTTTCCTG NNNNNNNNN`

Documentation
-------------
RNA library-level metadata
* [10X metadata](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-metrics)
* [Optimus metadata](https://broadinstitute.github.io/warp/docs/Pipelines/Optimus_Pipeline/Library-metrics)

RNA read-level metadata
* [10X molecule_info.h5](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-molecule-info)
* [Optimus .h5ad](https://broadinstitute.github.io/warp/docs/Pipelines/Optimus_Pipeline/Loom_schema)


