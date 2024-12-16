recon-count.jl
--------------
FASTQ parsing:
* `reads`: total number of reads
* `R1/R2_tooshort`: reads where R1/R2 length is shorter than the bead sequence structure
* `R1/R2_no_UP`: reads where the R1/R2 UP site does not match the expected sequence
* `R1/R2_GG_UP`: reads where the R1/R2 UP site is mostly "G"
* `R1/R2_N_UMI`: reads where the R1/R2 bead UMI contains an N
* `R1/R2_homopolymer_UMI`: reads where the R1/R2 bead UMI is highly degenerate
* `R1/R2_N_SB`: reads where the R1/R2 bead barcode contains an N
* `R1/R2_homopolymer_SB`: reads where the R1/R2 bead barcode is largely a homopolymer
* `reads/umis_filtered`: number of reads that pass the above filters - contain a clean (sb1, umi1, sb2, umi2)

Bead calling:
* `R1/R2_umicutoff`: number of umis a sb1/sb2 bead barcode needs to be called
* `R1/R2_barcodes`: number of called sb1/sb2 bead barcodes
* `R1/R2_exact`: number of umis belonging to a called sb1/sb2 bead barcode
* `umis_exact`: number of umis where both sb1 and sb2 exact match a called bead

Chimerism filtering:
* `R1/R2_chimeric`: number of umis removed for being chimeric in sb1/sb2
* `umis_chimeric`: number of umis where neither end is chimeric

Connection filter:
* `R1/R2_cxnfilter_cutoff`: number of connections above which a bead is filtered
* `R1/R2_cxnfilter_beads`: number of sb1/sb2 beads removed by the connection filter
* `R1/R2_cxnfilter`: number of umis containing a connection-filtered sb1/sb2 bead
* `umis_cxnfilter`: total number of umis removed by the connection filter

Final output results:
* `umis_final`: number of umis remaining after all matching+filtering steps
* `connections_final`: number of connections remaining in the final object

knn.jl
------

recon.py
--------
