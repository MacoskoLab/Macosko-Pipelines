Overview
--------
<pre>
FASTQs                   Diffusion matrix          KNN matrix             (sb,x,y)
R1s/R2s &rarr; recon-count.jl &rarr; matrix.csv.gz &rarr; knn.py &rarr; knn2.npz &rarr; recon.py &rarr; Puck.csv  
</pre>

recon-count.jl
--------------

General steps:
* extract clean (sb1, umi1, sb2, umi2) tuples from the FASTQs
* call bead barcodes with read counts above the elbow plot inflection point
* remove chimeras
* aggregatve reads into UMIs
* filter out high-connection beads
* save output

**matrix.csv.gz**: the diffusion matrix
* 3-column format `(sb1_index, sb2_index, umi)`
* `sb1.txt.gz`, `sb2.txt.gz` contain the sb1, sb2 barcode sequences (1-indexed)

**QC.pdf**: a visual summary of the input/output data

**metadata.csv**: contains processing metrics (see below)

* FASTQ parsing:
    * `reads`: total number of reads
    * `R1/R2_tooshort`: reads where R1/R2 length is shorter than the bead sequence structure
    * `R1/R2_no_UP`: reads where the R1/R2 UP site does not match the expected sequence
    * `R1/R2_GG_UP`: reads where the R1/R2 UP site is mostly "G"
    * `R1/R2_N_UMI`: reads where the R1/R2 bead UMI contains an N
    * `R1/R2_homopolymer_UMI`: reads where the R1/R2 bead UMI is highly degenerate
    * `R1/R2_N_SB`: reads where the R1/R2 bead barcode contains an N
    * `R1/R2_homopolymer_SB`: reads where the R1/R2 bead barcode is largely a homopolymer
    * `reads/umis_filtered`: number of reads/umis that pass the above filters
* Bead calling:
    * `R1/R2_umicutoff`: number of umis a sb1/sb2 bead barcode needs to be called
    * `R1/R2_barcodes`: number of called sb1/sb2 bead barcodes
    * `R1/R2_exact`: number of umis belonging to a called sb1/sb2 bead barcode
    * `umis_exact`: number of umis where both sb1 and sb2 exact match a called bead
* Chimerism filtering:
    * `R1/R2_chimeric`: number of umis removed for being chimeric in sb1/sb2
    * `umis_chimeric`: number of umis where either end is chimeric
* Connection filter:
    * `R1/R2_cxnfilter_z`: the z-score used to filter high-connection beads
    * `R1/R2_cxnfilter_cutoff`: number of connections above which a bead is filtered
    * `R1/R2_cxnfilter_beads`: number of sb1/sb2 beads removed by the connection filter
    * `R1/R2_cxnfilter`: number of umis containing a connection-filtered sb1/sb2 bead
    * `umis_cxnfilter`: total number of umis removed by the connection filter
* Final output results:
    * `umis_final`: number of umis remaining after all matching+filtering steps
    * `connections_final`: number of connections remaining in the final object
* General info:
    * `R1/R2_beadtype`: V10/V17 for R1, V15/V16 for R2
    * `downsampling_pct`: the % of reads retained from the original FASTQs

knn.py
------

Input:
* `in_dir`: path to dir with matrix.csv.gz
* `bead` (default 2)
* `n_neighbors` (default 150)

Output:
* sparse KNN matrix (`knn1.npz` or `knn2.npz`) in cosine distance

recon.py
--------
TODO
