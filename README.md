# Workflows

## bcl2fastq

### Inputs
bcl: gs:// path to BCL  
samplesheet: gs:// path to samplesheet  
technique: "cellranger" or "cellranger-arc" (TODO bcl2fastq)

### Outputs
/fastqs: output fastqs  
/logs: mkfastq.log and mkfastq.usage

# Docker Images
## "us-central1-docker.pkg.dev/velina-208320/docker-bcl2fastq/img:latest"
bcl2fastq2 v2.20  
Cell Ranger 8.0.0

# Test BCLs
Normal: 230521_SL-NXA_2112_AH3FLKBGXT (4 lanes)
