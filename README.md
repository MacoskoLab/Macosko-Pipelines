# Workflows

## bcl2fastq

**Inputs**  
bcl: gs:// path to BCL  
samplesheet: gs:// path to samplesheet  
technique: "cellranger" or "cellranger-arc" (TODO bcl2fastq)

**Commands**  
cellranger mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui  
cellranger-arc mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui  
rm -rf mkfastq/MAKE_FASTQS_CS  

**Outputs**  
/fastqs: output fastqs  
/logs: mkfastq.log and mkfastq.usage

**Notes**
* memory: "64 GB", cpu: 8, disks: "local-disk {max(BCL*3,64)} HDD"  
* throws an error if the disk is >6TB (edit bcl2fastq.wdl to increase cap)
* throws an error if the fastq directory already exists

## cellranger-count

**Inputs**  
bcl: gs:// path to BCL  
samplesheet: gs:// path to samplesheet  
technique: "cellranger" or "cellranger-arc" (TODO bcl2fastq)

**Commands**  
cellranger mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui  
cellranger-arc mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui  
rm -rf mkfastq/MAKE_FASTQS_CS  

**Outputs**  
/fastqs: output fastqs  
/logs: mkfastq.log and mkfastq.usage

**Notes**
* memory: "64 GB", cpu: 8, disks: "local-disk {max(BCL*3,64)} HDD"  
* throws an error if the disk is >6TB (edit bcl2fastq.wdl to increase cap)
* throws an error if the fastq directory already exists

# Docker Images
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
