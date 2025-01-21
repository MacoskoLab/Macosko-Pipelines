Seurat Metadata (TODO)
---------------
* orig.ident: `Slide-tags`
* nCount_RNA: total number of molecules detected within a cell
* nFeature_RNA: number of genes detected in each cell
* logumi: `log10(nCount_RNA+1)`
* cb: cell barcode, with trailing `-[0-9]*` removed
* cb_index: 1-indexed cell barcode enumeration
* percent.mt: `PercentageFeatureSet` with `pattern="^(MT-|mt-)"`
* pct.intronic:
* RNA_snn_res.0.8: `FindNeighbors` with `dims=1:30`, then `FindClusters` with `resolution=0.8`
* seurat_clusters: same as above

Seurat Misc
-----------
RNA_metrics?





