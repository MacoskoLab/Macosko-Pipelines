version 1.0

## Demultiplexing Genetic Samples Pipeline
##
## This workflow processes pooled scRNA-seq data to perform demultiplexing with functionality of removing mislabels/low depth samples using genetic data.
## It uses cellsnp-lite for Genotyping Bi-Allelic SNPs on Single Cells, bcftools for SNP or sample subsetting, and vireo for demultiplexing pooled scRNA-seq data with genotype reference.
##
## Requirements/expectations:
## - BAM file from each CellRanger/bender run.
## - Barcode file for cells.
## - QC'ed VCF File (plink2 --maf 0.2 --mac 20  --chr 1-22 --max-alleles 2 --vcf-min-gq 30 --vcf-min-dp 20 --export vcf 'bgz' 'vcf-dosage=HDS-force')
##
## Outputs:
## - Vireo output directories for demultiplexed libraries.
##
workflow DemultiplexSamples_removemis {
  input {
    Array[String] LIBNAMES
    String gsbucket_cellbender = "gs://fc-secure-e9137293-6b1d-47ad-ad37-fc170e7bf33d/libraries"
    String gsbucket_GVCF="gs://fc-secure-e9137293-6b1d-47ad-ad37-fc170e7bf33d/Stevens_Macosko_iNPH_WGS_filtered.vcf.gz"
    String gsbucket_GVCF_index="gs://fc-secure-e9137293-6b1d-47ad-ad37-fc170e7bf33d/Stevens_Macosko_iNPH_WGS_filtered.vcf.gz.tbi"
    String demx_docker = "us-central1-docker.pkg.dev/velina-208320/macosko-vireo/img"
    String thisBucket="gs://fc-secure-e9137293-6b1d-47ad-ad37-fc170e7bf33d"
    String thisBucket_libfolder= "~{thisBucket}/RESULTS"
    # Needed for if statements below, just ignore
    Boolean do_cellsnp = true
    Boolean do_subsetvcf_pileup  = true
    Boolean do_vireo_nodoublet    = true
    Boolean do_subsetvcf_samples     = true
    Boolean do_vireo_withdoublet = true



    File? noneFile
  }

  scatter (thisLibName in LIBNAMES) {
    
      if(do_cellsnp){
      call Run_cellsnp {
          input:
             thisLibName = thisLibName,
             gsbucket_cellbender = gsbucket_cellbender,
             gsbucket_GVCF = gsbucket_GVCF,
             gsbucket_GVCF_index=gsbucket_GVCF_index,
             thisBucket=thisBucket,
             thisBucket_libfolder=thisBucket_libfolder,
             docker = demx_docker
      }
      }
      if(do_subsetvcf_pileup){
      call Run_subsetvcf_pileup {
          input: 
        
             thisLibName = thisLibName,
             gsbucket_GVCF = gsbucket_GVCF,
             gsbucket_GVCF_index = gsbucket_GVCF_index,
             thisBucket=thisBucket,
             thisBucket_libfolder=thisBucket_libfolder,
             forceOrder = if do_cellsnp then Run_cellsnp.done else noneFile,
             docker = demx_docker
      }
      }
      if(do_vireo_nodoublet){
   call Run_vireo_nodoublet {
          input:
          
             thisLibName = thisLibName,
             thisBucket=thisBucket,
             thisBucket_libfolder=thisBucket_libfolder,
             forceOrder = if do_subsetvcf_pileup then Run_subsetvcf_pileup.done else noneFile,
             docker = demx_docker
      }
      }
      if(do_subsetvcf_samples){
   call Run_subsetvcf_samples {
          input: 
             thisLibName = thisLibName,
             gsbucket_GVCF = gsbucket_GVCF,
             gsbucket_GVCF_index=gsbucket_GVCF_index,
             thisBucket=thisBucket,
             thisBucket_libfolder=thisBucket_libfolder,
             forceOrder = if do_vireo_nodoublet then Run_vireo_nodoublet.done else noneFile,
             docker = demx_docker
      }
      }
      if(do_vireo_withdoublet){
   call Run_vireo_withdoublet {
          input:
             thisLibName = thisLibName,
             thisBucket=thisBucket,
             thisBucket_libfolder=thisBucket_libfolder,
             forceOrder = if do_subsetvcf_samples then Run_subsetvcf_samples.done else noneFile,
             docker = demx_docker
      }
      }
  }
}



task Run_cellsnp {
  input {
    String thisLibName
    String gsbucket_cellbender
    String gsbucket_GVCF
    String gsbucket_GVCF_index
    String thisBucket_libfolder
    String thisBucket
    String docker
    Int Memory = 16
    Int Disk = 120
    Int CPU = 4
    File? forceOrder
    Int? Preemtible = 1

    Float MakeSureBam = size("~{gsbucket_cellbender}/~{thisLibName}/outs/possorted_genome_bam.bam")
    Float MakeSureBamIdx = size("~{gsbucket_cellbender}/~{thisLibName}/outs/possorted_genome_bam.bam.bai")
    Float MakeSureFilterVCF = size("~{gsbucket_GVCF}")
    Float MakeSureBcds = size("~{gsbucket_cellbender}/~{thisLibName}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
    Float MakeSureFilterVCF_index = size("~{gsbucket_GVCF_index}")
  }

  command <<<
    touch /cromwell_root/usage.csv; dstat -t --cpu --mem --disk --io --freespace > /cromwell_root/usage.csv 2>&1 &

    cat > /tmp/RUNNER << 'EOF'
. /root/.bashrc; eval "$(/bin/micromamba shell hook -s posix)"; micromamba activate base;
mkdir torun; cd torun

gcloud storage cp ~{gsbucket_cellbender}/~{thisLibName}/outs/possorted_genome_bam.bam .
gcloud storage cp ~{gsbucket_cellbender}/~{thisLibName}/outs/possorted_genome_bam.bam.bai .
gcloud storage cp ~{gsbucket_GVCF} filtered.vcf.gz
gcloud storage cp ~{gsbucket_cellbender}/~{thisLibName}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz .
gcloud storage cp ~{gsbucket_GVCF_index} filtered.vcf.gz.tbi

mkdir cellsnp_out
echo "Running cellSNP-lite"
cellsnp-lite                       \
    -s possorted_genome_bam.bam    \
    -b barcodes.tsv.gz             \
    -R filtered.vcf.gz             \
    -O cellsnp_out                 \
    -p ~{CPU}
    # Check and compress cellSNP.base.vcf if .gz is missing
    if [ ! -s cellsnp_out/cellSNP.base.vcf.gz ]; then
      if [ -s cellsnp_out/cellSNP.base.vcf ]; then
        echo "Compressing cellSNP.base.vcf"
        bgzip -c cellsnp_out/cellSNP.base.vcf > cellsnp_out/cellSNP.base.vcf.gz
        cat /tmp/RUNNER > /cromwell_root/goodTask
      else
        echo "cellSNP.base.vcf.gz and cellSNP.base.vcf missing, cannot proceed."
        exit 1
      fi
    fi

export OUTPATH="~{thisBucket_libfolder}/~{thisLibName}/01_cellSNP" ;

date > .folder 
gcloud storage cp .folder ${OUTPATH}/.folder
gcloud storage cp cellsnp_out/cellSNP.tag.AD.mtx ${OUTPATH}/cellSNP.tag.AD.mtx
gcloud storage cp cellsnp_out/cellSNP.tag.DP.mtx ${OUTPATH}/cellSNP.tag.DP.mtx
gcloud storage cp cellsnp_out/cellSNP.tag.OTH.mtx ${OUTPATH}/cellSNP.tag.OTH.mtx
gcloud storage cp cellsnp_out/cellSNP.samples.tsv ${OUTPATH}/cellSNP.samples.tsv
gcloud storage cp cellsnp_out/cellSNP.base.vcf.gz ${OUTPATH}/cellSNP.base.vcf.gz
date
EOF

    bash /tmp/RUNNER |& ts

    touch /cromwell_root/usage.csv; pkill dstat; sleep 5s; pkill dstat; pkill dstat;
    true
  >>>
  output {
    File usage = "/cromwell_root/usage.csv"
    File done = "/cromwell_root/goodTask"
    }
  runtime {
    
    docker: docker
    memory: "~{Memory} GB"
    disks: "local-disk ~{Disk} HDD"
    cpu: "~{CPU}"
    preemtible: Preemtible
    zones: ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
  }
}





task Run_subsetvcf_pileup {
  input {
    String thisLibName
    String gsbucket_GVCF
    String gsbucket_GVCF_index
    String thisBucket_libfolder
    String thisBucket
    String docker
    Int Memory = 8
    Int Disk = 10
    Int CPU = 2
    File? forceOrder
    Int? Preemtible = 1

 Float MakeSureFilterVCF = size("~{gsbucket_GVCF}")

  }

  command <<<
    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:37.27.24.244:9008
    touch /cromwell_root/usage.csv; dstat -t --cpu --mem --disk --io --freespace > /cromwell_root/usage.csv 2>&1 &


    # Need single quote, was being weird with actually expanding and then failing
    cat > /tmp/RUNNER << 'EOF'

. /root/.bashrc; eval "$(/bin/micromamba shell hook -s posix)"; micromamba activate base;
mkdir torun; cd  torun

gcloud storage cp ~{thisBucket_libfolder}/~{thisLibName}/01_cellSNP/cellSNP.base.vcf.gz .
gcloud storage cp  ~{gsbucket_GVCF} filtered.vcf.gz 
gcloud storage cp ~{gsbucket_GVCF_index} filtered.vcf.gz.tbi


bcftools view  filtered.vcf.gz -R cellSNP.base.vcf.gz -Ov -o pileup_filtered.vcf.gz

if [ ! -s pileup_filtered.vcf.gz ]; then
    echo "Failed. VCF empty."
    exit 1
else
    echo "VCF subsetting done"
    cat /tmp/RUNNER > /cromwell_root/goodTask
fi


export OUTPATH="~{thisBucket_libfolder}/~{thisLibName}/02_VCFSUBSET1"



# Gcloud changes what it does depending if a folder already exists. So just force it with a blan file
date > .folder 
gcloud storage cp .folder ${OUTPATH}/.folder
gcloud storage cp pileup_filtered.vcf.gz  ${OUTPATH}/pileup_filtered.vcf.gz
date
EOF

    bash /tmp/RUNNER |& ts

    touch /cromwell_root/usage.csv; pkill dstat; sleep 5s; pkill dstat; pkill dstat;
    true

  >>>
  output {
    File usage = "/cromwell_root/usage.csv"
    # Helps prevent call caching, I think, if make sure this is only created on done
  	File done = "/cromwell_root/goodTask"

  }
  runtime {
   
    docker: docker
    memory: "~{Memory} GB"
    disks: "local-disk ~{Disk} HDD"
    cpu: "~{CPU}"
    preemtible: Preemtible
    zones: ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
  }
}


task Run_vireo_nodoublet {
  input {
    String thisLibName
    String thisBucket_libfolder
    String thisBucket
    String docker
    Int Memory = 80
    Int Disk = 40
    Int CPU = 4
    File? forceOrder
    Int? Preemtible = 1
    Float MakeSureFilterCELL = size("~{thisBucket_libfolder}/~{thisLibName}/01_cellSNP/cellSNP.base.vcf.gz")
    Float MakeSureFilterVCF = size("~{thisBucket_libfolder}/~{thisLibName}/02_VCFSUBSET1/pileup_filtered.vcf.gz")
  }

  command <<<
    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:37.27.24.244:9008
    touch /cromwell_root/usage.csv; dstat -t --cpu --mem --disk --io --freespace > /cromwell_root/usage.csv 2>&1 &


    # Need single quote, was being weird with actually expanding and then failing
    cat > /tmp/RUNNER << 'EOF'

. /root/.bashrc; eval "$(/bin/micromamba shell hook -s posix)"; micromamba activate base;
mkdir torun; cd  torun

gcloud storage cp -r ~{thisBucket_libfolder}/~{thisLibName}/01_cellSNP/ .
gcloud storage cp -r ~{thisBucket_libfolder}/~{thisLibName}/02_VCFSUBSET1/pileup_filtered.vcf.gz . 


mkdir vireo_all_out
vireo -c 01_cellSNP/ -o vireo_all_out -d pileup_filtered.vcf.gz -t GT --noDoublet 

if [ ! -s "vireo_all_out/summary.tsv" ]; then
    echo "Vireo without doublet detection failed."
    exit 1
else
    echo "Vireo without doublet detection is done"
    cat /tmp/RUNNER > /cromwell_root/goodTask
fi


cat vireo_all_out/summary.tsv |awk 'NR>1 { if($2>50) print $1}'|grep -v "unassigned" > samples.txt
export OUTPATH="~{thisBucket_libfolder}/~{thisLibName}/03_VIREO1"

# Gcloud changes what it does depending if a folder already exists. So just force it with a blan file
date > .folder 
gcloud storage cp .folder ${OUTPATH}/
gcloud storage cp -r vireo_all_out  ${OUTPATH}/
gcloud storage cp samples.txt ${OUTPATH}/

date
EOF

    bash /tmp/RUNNER |& ts

    touch /cromwell_root/usage.csv; pkill dstat; sleep 5s; pkill dstat; pkill dstat;
  true

  >>>
  output {
    File usage = "/cromwell_root/usage.csv"
    # Helps prevent call caching, I think, if make sure this is only created on done
  	File done = "/cromwell_root/goodTask"
  }
  runtime {
  
    docker: docker
    memory: "~{Memory} GB"
    disks: "local-disk ~{Disk} HDD"
    cpu: "~{CPU}"
    preemtible: Preemtible
    zones: ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
  }
}


task Run_subsetvcf_samples {
  input {
    String thisLibName
    String gsbucket_GVCF
    String gsbucket_GVCF_index
    String thisBucket_libfolder
    String thisBucket
    String docker
    Int Memory = 8
    Int Disk = 10
    Int CPU = 2
    File? forceOrder
    Int? Preemtible = 1

 Float MakeSureFilterVCF = size("~{gsbucket_GVCF}")
 Float MakeSureFilterVCF_index = size("~{gsbucket_GVCF_index}")
  }

  command <<<
    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:37.27.24.244:9008
    touch /cromwell_root/usage.csv; dstat -t --cpu --mem --disk --io --freespace > /cromwell_root/usage.csv 2>&1 &


    # Need single quote, was being weird with actually expanding and then failing
    cat > /tmp/RUNNER << 'EOF'

. /root/.bashrc; eval "$(/bin/micromamba shell hook -s posix)"; micromamba activate base;
mkdir torun; cd  torun

gcloud storage cp ~{thisBucket_libfolder}/~{thisLibName}/01_cellSNP/cellSNP.base.vcf.gz .
gcloud storage cp  ~{gsbucket_GVCF} filtered.vcf.gz 
gcloud storage cp  ~{gsbucket_GVCF_index} filtered.vcf.gz.tbi
gcloud storage cp ~{thisBucket_libfolder}/~{thisLibName}/03_VIREO1/samples.txt .


bcftools view  -S samples.txt filtered.vcf.gz -R cellSNP.base.vcf.gz  -Ov -o subsetted_filtered.vcf.gz
if [ ! -s "subsetted_filtered.vcf.gz" ]; then
    echo "vcftools subsetting vireo samples failed"
    exit 1
else
    echo "vcftools subsetting vireo samples done"
    cat /tmp/RUNNER > /cromwell_root/goodTask
fi

export OUTPATH="~{thisBucket_libfolder}/~{thisLibName}/04_VCFSUBSET2"

# Gcloud changes what it does depending if a folder already exists. So just force it with a blan file
date > .folder 
gcloud storage cp .folder ${OUTPATH}/.folder
gcloud storage cp subsetted_filtered.vcf.gz  ${OUTPATH}/subsetted_filtered.vcf.gz
date
EOF

    bash /tmp/RUNNER |& ts

    touch /cromwell_root/usage.csv; pkill dstat; sleep 5s; pkill dstat; pkill dstat;
   true

  >>>
  output {
    File usage = "/cromwell_root/usage.csv"
    # Helps prevent call caching, I think, if make sure this is only created on done
  	File done = "/cromwell_root/goodTask"
  }
  runtime {
    docker: docker
    memory: "~{Memory} GB"
    disks: "local-disk ~{Disk} HDD"
    cpu: "~{CPU}"
    preemtible: Preemtible
    zones: ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
  }
}


task Run_vireo_withdoublet {
  input {
    String thisLibName
    String thisBucket_libfolder
    String thisBucket
    String docker
    Int Memory = 80
    Int Disk = 40
    Int CPU = 4
    File? forceOrder
    Int? Preemtible = 1
   Float MakeSureFinalVCF = size("~{thisBucket_libfolder}/~{thisLibName}/04_VCFSUBSET2/subsetted_filtered.vcf.gz")
  }

  command <<<
    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:37.27.24.244:9008
    touch /cromwell_root/usage.csv; dstat -t --cpu --mem --disk --io --freespace > /cromwell_root/usage.csv 2>&1 &


    # Need single quote, was being weird with actually expanding and then failing
    cat > /tmp/RUNNER << 'EOF'

. /root/.bashrc; eval "$(/bin/micromamba shell hook -s posix)"; micromamba activate base;
mkdir torun; cd  torun

gcloud storage cp -r ~{thisBucket_libfolder}/~{thisLibName}/01_cellSNP/ .
gcloud storage cp -r ~{thisBucket_libfolder}/~{thisLibName}/04_VCFSUBSET2/subsetted_filtered.vcf.gz . 

mkdir vireo_final_out
vireo -c 01_cellSNP/ -o vireo_final_out -d subsetted_filtered.vcf.gz -t GT 

if [ ! -s "vireo_final_out/summary.tsv" ]; then
    echo "Vireo with doublet detection failed."
    exit 1
else
    echo "Vireo with doublet detection is done"
    cat /tmp/RUNNER > /cromwell_root/goodTask
fi


export OUTPATH="~{thisBucket_libfolder}/~{thisLibName}/05_VIREO"

# Gcloud changes what it does depending if a folder already exists. So just force it with a blan file
date > .folder 
gcloud storage cp .folder ${OUTPATH}/
gcloud storage cp -r vireo_final_out  ${OUTPATH}/




date
EOF

    bash /tmp/RUNNER |& ts

    touch /cromwell_root/usage.csv; pkill dstat; sleep 5s; pkill dstat; pkill dstat;
 true
  >>>
  output {
    File usage = "/cromwell_root/usage.csv"
    # Helps prevent call caching, I think, if make sure this is only created on done
  	File done = "/cromwell_root/goodTask"
  }
  runtime {
    
    docker: docker
    memory: "~{Memory} GB"
    disks: "local-disk ~{Disk} HDD"
    cpu: "~{CPU}"
    preemtible: Preemtible
    zones: ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
  }
}
