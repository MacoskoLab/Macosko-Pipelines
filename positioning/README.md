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

RNA source

call-CellMetrics/attempt-2/Test_Optimus_gex.cell-metrics.csv.gz

call-GeneMetrics/Test_Optimus_gex.gene-metrics.csv.gz

call-ReferenceCheck/cacheCopy/reference_version.txt

call-OptimusH5adGeneration/Test_Optimus_gex.h5ad
call-OptimusH5adGeneration/Test_Optimus_gex__library_metrics.csv

call-MergeBam/Test_Optimus_gex.bam

call-MergeStarOutputs/Test_Optimus_gex.star_metrics.tar
call-MergeStarOutputs/Test_Optimus_gex_sparse_counts_row_index.npy
call-MergeStarOutputs/Test_Optimus_gex_sparse_counts_col_index.npy
call-MergeStarOutputs/Test_Optimus_gex_filtered_mtx_files.tar
call-MergeStarOutputs/Test_Optimus_gex.mtx_files.tar
