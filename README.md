# Workflows

## bcl2fastq

#### Inputs
bcl: gs:// path to BCL  
samplesheet: gs:// path to samplesheet  
technique: "cellranger" or "cellranger-arc" (TODO bcl2fastq)

#### Commands
cellranger mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui  
cellranger-arc mkfastq --run=BCL --id=mkfastq --csv=Indexes.csv --disable-ui  
rm -rf mkfastq/MAKE_FASTQS_CS  

#### Outputs
/fastqs: output fastqs  
/logs: mkfastq.log and mkfastq.usage

#### Notes
memory: "64 GB", cpu: 8, disks: "local-disk {max(BCL*3,64)} LOCAL"


# Docker Images
**us-central1-docker.pkg.dev/velina-208320/docker-bcl2fastq/img:latest**  
bcl2fastq2 v2.20  
Cell Ranger 8.0.0



# Test BCLs
Normal: 230521_SL-NXA_2112_AH3FLKBGXT (4 lanes, 24.06 GiB)
