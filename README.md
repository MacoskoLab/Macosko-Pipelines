Workflows
=========

bcl2fastq
---------

**Inputs**  
* bcl: gs:// path to BCL  
* samplesheet: gs:// path to samplesheet  
* technique: "cellranger" or "cellranger-arc" or "cellranger-atac" or "bcl2fastq"
* fastq_output_path (optional): gs:// path to write fastqs (default fastqs/{basename(bcl)})  
* log_output_path (optional): gs:// path to write logs (default logs/{basename(bcl)})

**Commands**  
* cellranger mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui  
* cellranger-arc mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui  
* cellranger-atac mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui
* bcl2fastq --runfolder-dir BCL --input-dir BCL/Data/Intensities/BaseCalls --output-dir mkfastq --sample-sheet Indexes.csv --create-fastq-for-index-reads  
* rm -rf mkfastq/MAKE_FASTQS_CS  

**Outputs**  
* /fastqs: output fastqs  
* /logs: mkfastq.log and mkfastq.usage

**Notes**
* memory: "64 GB", cpu: 8, disks: "local-disk {max(BCL*3,128)} SSD"  
* throws an error if the disk is >6TB (edit bcl2fastq.wdl to increase cap)  
* throws an error if the fastq directory already exists

cellranger-count
----------------

**Inputs**  
* fastq_path: gs:// path to the fastq folder  
* sample: fastq filename prefix to select (specified in the sample sheet supplied to the FASTQ generation software)  
* reference: gs:// path to the transcriptome  
* technique: "cellranger" or "cellranger-atac"  
* lanes (optional): Array[Int] of lanes to subset (default is [], meaning all lanes)  
* count_output_path (optional): gs:// path to write outs (default cellranger-count/{basename(fastq_path)})  
* log_output_path (optional): gs:// path to write logs (default logs/{basename(fastq_path)})

**Commands**  
* cellranger count --id={id} --transcriptome=reference --fastqs=fastqs --sample={sample} --create-bam=true --include-introns=true --nosecondary --disable-ui  
* cellranger-atac count --id={id} --reference=reference --fastqs=fastqs --disable-ui  
* rm -rf {id}/SC_RNA_COUNTER_CS

**Outputs**  
* /cellranger-count: output cellranger results  
* /logs: count-{id}.out, count-{id}.err, count-{id}.usage  

**Notes**
* memory: "64 GB", cpu: 8, disks: "local-disk {max(fastqs\*6+20,128)} SSD"  
* throws an error if the disk is >6TB (edit cellranger-count.wdl to increase cap)
* throws an error if the counts directory already exists
* cellranger expects the fastqs to be named as [sample]\_S[number]\_L00[lane]\_[R1/R2/I1/I2]\_001.fastq.gz
* the output folder is {id}: equal to {sample} if all lanes are used, {sample}_{lanes} otherwise

spatial-count
----------------

**Inputs**  
* fastq_path: gs:// path to the fastq folder  
* sample: fastq filename prefix to select (specified in the sample sheet supplied to the FASTQ generation software)  
* pucks: array of gs:// paths to each puck  
* lanes (optional): Array[Int] of lanes to subset (default is [], meaning all lanes)  
* count_output_path (optional): gs:// path to write outs (default spatial-count/{basename(fastq_path)})  
* log_output_path (optional): gs:// path to write logs (default logs/{basename(fastq_path)})

**Commands**  
* julia spatial-count.jl fastqs pucks  

**Outputs**  
* /spatial-count: SBcounts.h5  
* /logs: count-{id}.out, count-{id}.err, count-{id}.usage  

**Notes**
* memory: "{max(fastqs\*2.5,64)} GB", cpu: 1, disks: "local-disk {max(fastqs\*2.5,64)} SSD"  
* throws an error if the disk is >256GB (edit spatial-count.wdl to increase cap)
* WDL expects the fastqs to be named as [sample]\_S[number]\_L00[lane]\_[R1/R2/I1/I2]\_001.fastq.gz
* the output folder is {id}: equal to {sample} if all lanes are used, {sample}_{lanes} otherwise

Input→Output
====================
* bcl2fastq: bcls, samplesheets → fastqs  
* cellranger-count: fastqs, references → cellranger-count  
* spatial-count: fastqs, pucks → spatial-count
* reconstruction: fastqs → reconstruction  

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
Python 2.7.16  
Python3 3.7.3  
Julia 1.8.5  
R 4.3.3  

**us-central1-docker.pkg.dev/velina-208320/docker-recon/img:latest**  
Python 3.11  
CUDA 12.2  

# Test BCLs
Normal: 230521_SL-NXA_2112_AH3FLKBGXT (4 lanes, 24.06 GiB)  
Multiome GEX: 230824_VH01286_38_AACVHVJM5 (1 lane, 23G)  
Multiome ATAC: 230827_VH01286_40_AACVHVLM5 (1 lane, 20G)  
Slide-Tags: 240329_SL-NTA_74_AAFJ3LYM5 (1 lane, 42.69 GiB)

Links
=====
Cell Ranger downloads:  
https://www.10xgenomics.com/support/software/cell-ranger/downloads  
Cell Ranger ATAC downloads:  
https://support.10xgenomics.com/single-cell-atac/software/downloads/latest  
Cell Ranger ARC downloads:  
https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/downloads/latest  

Helpful commands
================
.h5 commands
------------
* h5dump -n SBcounts.h5: list all contents
* h5dump -d /metadata/num_reads SBcounts.h5: print a specific dataset
