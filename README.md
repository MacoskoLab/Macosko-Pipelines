Workflows
=========

bcl2fastq
---------

**Inputs**  
bcl: gs:// path to BCL  
samplesheet: gs:// path to samplesheet  
technique: "cellranger" or "cellranger-arc" or "cellranger-atac" (TODO bcl2fastq)  
fastq_output_path (optional): gs:// path to write fastqs (default fastqs/{basename(bcl)})  
log_output_path (optional): gs:// path to write logs (default logs/{basename(bcl)})

**Commands**  
cellranger mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui  
cellranger-arc mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui  
cellranger-atac mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui  
rm -rf mkfastq/MAKE_FASTQS_CS  

**Outputs**  
/fastqs: output fastqs  
/logs: mkfastq.log and mkfastq.usage

**Notes**
* memory: "64 GB", cpu: 8, disks: "local-disk {max(BCL*3,128)} SSD"  
* throws an error if the disk is >6TB (edit bcl2fastq.wdl to increase cap)  
* throws an error if the fastq directory already exists

cellranger-count
----------------

**Inputs**  
fastqs: gs:// path to the fastq folder  
sample: fastq filename prefix to select (specified in the sample sheet supplied to the FASTQ generation software)  
reference: gs:// path to the transcriptome  
technique: "cellranger" or "cellranger-atac"  
lanes (optional): Array[Int] of lanes to subset (default is [], meaning all lanes)  
count_output_path (optional): gs:// path to write outs (default cellranger-count/{basename(fastqs)})  
log_output_path (optional): gs:// path to write logs (default logs/{basename(fastqs)})

**Commands**  
cellranger count  
cellranger-arc count  
rm -rf {id}/SC_RNA_COUNTER_CS

**Outputs**  
/fastqs: output fastqs  
/logs: mkfastq.log and mkfastq.usage

**Notes**
* memory: "64 GB", cpu: 8, disks: "local-disk {max(fastqs*6+20,128)} SSD"  
* throws an error if the disk is >6TB (edit count.wdl to increase cap)
* throws an error if the counts directory already exists
* cellranger expects the fastqs to be named as [sample]\_S[number]\_L00[lane]\_[R1/R2/I1/I2]\_001.fastq.gz
* the output folder is {id}: equal to {sample} if all lanes are used, {sample}_{lanes} otherwise

Input→Output
====================
bcl2fastq: bcls, samplesheets → fastqs  
cellranger-count: fastqs, references → cellranger-count  

Docker Images
=============
**us-central1-docker.pkg.dev/velina-208320/docker-bcl2fastq/img:latest**  
bcl2fastq2 v2.20  
Cell Ranger 8.0.0  
Cell Ranger ATAC 2.1.0  
Cell Ranger ARC 2.0.2  

**us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest**  
Cell Ranger 8.0.0  
Cell Ranger ATAC 2.1.0  
Cell Ranger ARC 2.0.2  
Python3 3.7.3  
Julia 1.8.5  
R 4.3.3

# Test BCLs
Normal: 230521_SL-NXA_2112_AH3FLKBGXT (4 lanes, 24.06 GiB)  
Multiome GEX: 230824_VH01286_38_AACVHVJM5 (1 lane, 23G)  
Multiome ATAC: 230827_VH01286_40_AACVHVLM5 (1 lane, 20G)  

Links
=====
Cell Ranger downloads:  
https://www.10xgenomics.com/support/software/cell-ranger/downloads  
Cell Ranger ATAC downloads:  
https://support.10xgenomics.com/single-cell-atac/software/downloads/latest  
Cell Ranger ARC downloads:  
https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/downloads/latest  
