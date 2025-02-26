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
RNA_metadata
* [10X metadata](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-metrics)
* [Optimus metadata](https://broadinstitute.github.io/warp/docs/Pipelines/Optimus_Pipeline/Library-metrics)

Other documentation
* [10X molecule_info.h5](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-molecule-info)
* [Optimus .h5ad](https://broadinstitute.github.io/warp/docs/Pipelines/Optimus_Pipeline/Loom_schema)

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
    * R1_tooshort - number of reads with fewer than the allowable 
    * R2_tooshort -
    * UMI - filtering statistics on the UMI sequence
    * UP - filtering statistics on the UP
    * SB - filtering statistics on the spatial barcode sequence
    * SB_HD - location of spatial barcode fuzzy match

Bead structures
---------------
* `JJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJNNNNNNNVV`   (V10)
* `JJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJJNNNNNNNNN` (V17)
* `JJJJJJJJJJJJJJJCTGTTTCCTGNNNNNNNNN`           (V15)
* `JJJJJJJJJJJJJJJJJCTGTTTCCTGNNNNNNNNN`         (V16)

