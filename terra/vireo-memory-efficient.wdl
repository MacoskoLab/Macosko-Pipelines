version 1.0

workflow DemultiplexSamples_removemis {
  input {
    Array[String] LIBNAMES
    String gsbucket_cellbender = "gs://fc-secure-e9137293-6b1d-47ad-ad37-fc170e7bf33d/libraries"
    String gsbucket_GVCF       = "gs://fc-secure-e9137293-6b1d-47ad-ad37-fc170e7bf33d/Stevens_Macosko_iNPH_WGS_filtered.vcf.gz"
    String gsbucket_GVCF_index = "gs://fc-secure-e9137293-6b1d-47ad-ad37-fc170e7bf33d/Stevens_Macosko_iNPH_WGS_filtered.vcf.gz.tbi"
    String demx_docker         = "us-central1-docker.pkg.dev/velina-208320/macosko-vireo/img"
    String thisBucket          = "gs://fc-secure-e9137293-6b1d-47ad-ad37-fc170e7bf33d"
    String thisBucket_libfolder = "~{thisBucket}/RESULTS"

    Boolean do_cellsnp            = true
    Boolean do_subsetvcf_pileup   = true
    Boolean do_vireo_nodoublet    = true
    Boolean do_subsetvcf_samples  = true
    Boolean do_vireo_withdoublet  = true

    File? noneFile
  }

  scatter (thisLibName in LIBNAMES) {

    if (do_cellsnp) {
      call Run_cellsnp {
        input:
          thisLibName           = thisLibName,
          gsbucket_cellbender   = gsbucket_cellbender,
          gsbucket_GVCF         = gsbucket_GVCF,
          gsbucket_GVCF_index   = gsbucket_GVCF_index,
          thisBucket            = thisBucket,
          thisBucket_libfolder  = thisBucket_libfolder,
          docker                = demx_docker
      }
    }

    if (do_subsetvcf_pileup) {
      call Run_subsetvcf_pileup {
        input:
          thisLibName          = thisLibName,
          gsbucket_GVCF        = gsbucket_GVCF,
          gsbucket_GVCF_index  = gsbucket_GVCF_index,
          thisBucket           = thisBucket,
          thisBucket_libfolder = thisBucket_libfolder,
          docker               = demx_docker,
          forceOrder           = select_first([Run_cellsnp.done, noneFile])
      }
    }

    if (do_vireo_nodoublet) {
      call Run_vireo_nodoublet {
        input:
          thisLibName          = thisLibName,
          thisBucket           = thisBucket,
          thisBucket_libfolder = thisBucket_libfolder,
          docker               = demx_docker,
          forceOrder           = select_first([Run_subsetvcf_pileup.done, noneFile])
      }
    }

    if (do_subsetvcf_samples) {
      call Run_subsetvcf_samples {
        input:
          thisLibName          = thisLibName,
          gsbucket_GVCF        = gsbucket_GVCF,
          gsbucket_GVCF_index  = gsbucket_GVCF_index,
          thisBucket           = thisBucket,
          thisBucket_libfolder = thisBucket_libfolder,
          docker               = demx_docker,
          forceOrder           = select_first([Run_vireo_nodoublet.done, noneFile])
      }
    }

    if (do_vireo_withdoublet) {
      call Run_vireo_withdoublet {
        input:
          thisLibName          = thisLibName,
          thisBucket           = thisBucket,
          thisBucket_libfolder = thisBucket_libfolder,
          docker               = demx_docker,
          forceOrder           = select_first([Run_subsetvcf_samples.done, noneFile])
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

    Float MakeSureBam              = size("~{gsbucket_cellbender}/~{thisLibName}/outs/possorted_genome_bam.bam")
    Float MakeSureBamIdx           = size("~{gsbucket_cellbender}/~{thisLibName}/outs/possorted_genome_bam.bam.bai")
    Float MakeSureFilterVCF        = size("~{gsbucket_GVCF}")
    Float MakeSureBcds             = size("~{gsbucket_cellbender}/~{thisLibName}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
    Float MakeSureFilterVCF_index  = size("~{gsbucket_GVCF_index}")
  }

  command <<<
    # metrics log in working dir
    touch usage.csv
    if command -v dstat >/dev/null 2>&1; then
      dstat -t --cpu --mem --disk --io --freespace > usage.csv 2>&1 & DSTAT_PID=$!
    fi

    cat > /tmp/RUNNER << 'EOF'
set -euo pipefail
# Initialize micromamba without sourcing ~/.bashrc (avoids PS1 unbound var under set -u)
eval "$(/bin/micromamba shell hook -s posix)"
micromamba activate base

mkdir -p torun
cd torun

gcloud storage cp ~{gsbucket_cellbender}/~{thisLibName}/outs/possorted_genome_bam.bam .
gcloud storage cp ~{gsbucket_cellbender}/~{thisLibName}/outs/possorted_genome_bam.bam.bai .
gcloud storage cp ~{gsbucket_GVCF} filtered.vcf.gz
gcloud storage cp ~{gsbucket_cellbender}/~{thisLibName}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz .
gcloud storage cp ~{gsbucket_GVCF_index} filtered.vcf.gz.tbi

mkdir -p cellsnp_out
echo "Running cellSNP-lite"
cellsnp-lite \
  -s possorted_genome_bam.bam \
  -b barcodes.tsv.gz \
  -R filtered.vcf.gz \
  -O cellsnp_out \
  -p ~{CPU}

# ensure gz VCF exists; compress if needed and (optionally) index
if [ ! -s cellsnp_out/cellSNP.base.vcf.gz ]; then
  if [ -s cellsnp_out/cellSNP.base.vcf ]; then
    echo "Compressing cellSNP.base.vcf"
    bgzip -c cellsnp_out/cellSNP.base.vcf > cellsnp_out/cellSNP.base.vcf.gz
  else
    echo "cellSNP.base.vcf(.gz) missing, cannot proceed."
    exit 1
  fi
fi

if command -v bcftools >/dev/null 2>&1; then
  bcftools index -t cellsnp_out/cellSNP.base.vcf.gz || true
fi

# success marker one dir up (task working dir)
echo ok > ../goodTask
ls -l ..

export OUTPATH="~{thisBucket_libfolder}/~{thisLibName}/01_cellSNP"
date > .folder
gcloud storage cp .folder "${OUTPATH}/.folder"
gcloud storage cp cellsnp_out/cellSNP.tag.AD.mtx   "${OUTPATH}/cellSNP.tag.AD.mtx"
gcloud storage cp cellsnp_out/cellSNP.tag.DP.mtx   "${OUTPATH}/cellSNP.tag.DP.mtx"
gcloud storage cp cellsnp_out/cellSNP.tag.OTH.mtx  "${OUTPATH}/cellSNP.tag.OTH.mtx"
gcloud storage cp cellsnp_out/cellSNP.samples.tsv  "${OUTPATH}/cellSNP.samples.tsv"
gcloud storage cp cellsnp_out/cellSNP.base.vcf.gz  "${OUTPATH}/cellSNP.base.vcf.gz"
if [ -f cellsnp_out/cellSNP.base.vcf.gz.tbi ]; then
  gcloud storage cp cellsnp_out/cellSNP.base.vcf.gz.tbi "${OUTPATH}/cellSNP.base.vcf.gz.tbi"
fi
date
EOF

    if command -v ts >/dev/null 2>&1; then
      bash /tmp/RUNNER |& ts
    else
      bash /tmp/RUNNER
    fi

    if [ -n "${DSTAT_PID:-}" ]; then kill "$DSTAT_PID" 2>/dev/null || true; fi
  >>>

  output {
    File usage = "usage.csv"
    File done  = "goodTask"
  }

  runtime {
    docker: docker
    memory: "~{Memory} GB"
    disks: "local-disk ~{Disk} HDD"
    cpu: "~{CPU}"
    preemptible: Preemtible
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
    touch usage.csv
    if command -v dstat >/dev/null 2>&1; then
      dstat -t --cpu --mem --disk --io --freespace > usage.csv 2>&1 & DSTAT_PID=$!
    fi

    cat > /tmp/RUNNER << 'EOF'
set -euo pipefail
eval "$(/bin/micromamba shell hook -s posix)"
micromamba activate base

mkdir -p torun
cd torun

gcloud storage cp "~{thisBucket_libfolder}/~{thisLibName}/01_cellSNP/cellSNP.base.vcf.gz" .
gcloud storage cp "~{gsbucket_GVCF}" filtered.vcf.gz
gcloud storage cp "~{gsbucket_GVCF_index}" filtered.vcf.gz.tbi

# Compressed output + index
bcftools view filtered.vcf.gz -R cellSNP.base.vcf.gz -Oz -o pileup_filtered.vcf.gz
bcftools index -t pileup_filtered.vcf.gz

if [ ! -s pileup_filtered.vcf.gz ]; then
  echo "Failed. VCF empty."
  exit 1
else
  echo ok > ../goodTask
fi
ls -l ..

export OUTPATH="~{thisBucket_libfolder}/~{thisLibName}/02_VCFSUBSET1"
date > .folder
gcloud storage cp .folder "${OUTPATH}/.folder"
gcloud storage cp pileup_filtered.vcf.gz     "${OUTPATH}/pileup_filtered.vcf.gz"
gcloud storage cp pileup_filtered.vcf.gz.tbi "${OUTPATH}/pileup_filtered.vcf.gz.tbi"
date
EOF

    if command -v ts >/dev/null 2>&1; then
      bash /tmp/RUNNER |& ts
    else
      bash /tmp/RUNNER
    fi

    if [ -n "${DSTAT_PID:-}" ]; then kill "$DSTAT_PID" 2>/dev/null || true; fi
  >>>

  output {
    File usage = "usage.csv"
    File done  = "goodTask"
  }

  runtime {
    docker: docker
    memory: "~{Memory} GB"
    disks: "local-disk ~{Disk} HDD"
    cpu: "~{CPU}"
    preemptible: Preemtible
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
    Float MakeSureFilterVCF  = size("~{thisBucket_libfolder}/~{thisLibName}/02_VCFSUBSET1/pileup_filtered.vcf.gz")
  }

  command <<<
    touch usage.csv
    if command -v dstat >/dev/null 2>&1; then
      dstat -t --cpu --mem --disk --io --freespace > usage.csv 2>&1 & DSTAT_PID=$!
    fi

    cat > /tmp/RUNNER << 'EOF'
set -euo pipefail
eval "$(/bin/micromamba shell hook -s posix)"
micromamba activate base

mkdir -p torun
cd torun

gcloud storage cp -r "~{thisBucket_libfolder}/~{thisLibName}/01_cellSNP/" .
gcloud storage cp    "~{thisBucket_libfolder}/~{thisLibName}/02_VCFSUBSET1/pileup_filtered.vcf.gz" .
gcloud storage cp    "~{thisBucket_libfolder}/~{thisLibName}/02_VCFSUBSET1/pileup_filtered.vcf.gz.tbi" .

mkdir -p vireo_all_out
vireo -c 01_cellSNP/ -o vireo_all_out -d pileup_filtered.vcf.gz -t GT --noDoublet

if [ ! -s "vireo_all_out/summary.tsv" ]; then
  echo "Vireo without doublet detection failed."
  exit 1
else
  echo ok > ../goodTask
fi

awk 'NR>1 && $2>50 && $1!="unassigned" {print $1}' vireo_all_out/summary.tsv > samples.txt
ls -l ..

export OUTPATH="~{thisBucket_libfolder}/~{thisLibName}/03_VIREO1"
date > .folder
gcloud storage cp .folder "${OUTPATH}/.folder"
gcloud storage cp -r vireo_all_out "${OUTPATH}/"
gcloud storage cp samples.txt     "${OUTPATH}/"
date
EOF

    if command -v ts >/dev/null 2>&1; then
      bash /tmp/RUNNER |& ts
    else
      bash /tmp/RUNNER
    fi

    if [ -n "${DSTAT_PID:-}" ]; then kill "$DSTAT_PID" 2>/dev/null || true; fi
  >>>

  output {
    File usage = "usage.csv"
    File done  = "goodTask"
  }

  runtime {
    docker: docker
    memory: "~{Memory} GB"
    disks: "local-disk ~{Disk} HDD"
    cpu: "~{CPU}"
    preemptible: Preemtible
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

    Float MakeSureFilterVCF       = size("~{gsbucket_GVCF}")
    Float MakeSureFilterVCF_index = size("~{gsbucket_GVCF_index}")
  }

  command <<<
    touch usage.csv
    if command -v dstat >/dev/null 2>&1; then
      dstat -t --cpu --mem --disk --io --freespace > usage.csv 2>&1 & DSTAT_PID=$!
    fi

    cat > /tmp/RUNNER << 'EOF'
set -euo pipefail
eval "$(/bin/micromamba shell hook -s posix)"
micromamba activate base

mkdir -p torun
cd torun

gcloud storage cp "~{thisBucket_libfolder}/~{thisLibName}/01_cellSNP/cellSNP.base.vcf.gz" .
gcloud storage cp "~{gsbucket_GVCF}"       filtered.vcf.gz
gcloud storage cp "~{gsbucket_GVCF_index}" filtered.vcf.gz.tbi
gcloud storage cp "~{thisBucket_libfolder}/~{thisLibName}/03_VIREO1/samples.txt" .

bcftools view -S samples.txt filtered.vcf.gz -R cellSNP.base.vcf.gz -Oz -o subsetted_filtered.vcf.gz
bcftools index -t subsetted_filtered.vcf.gz

if [ ! -s "subsetted_filtered.vcf.gz" ]; then
  echo "bcftools subsetting vireo samples failed"
  exit 1
else
  echo ok > ../goodTask
fi
ls -l ..

export OUTPATH="~{thisBucket_libfolder}/~{thisLibName}/04_VCFSUBSET2"
date > .folder
gcloud storage cp .folder "${OUTPATH}/.folder"
gcloud storage cp subsetted_filtered.vcf.gz     "${OUTPATH}/subsetted_filtered.vcf.gz"
gcloud storage cp subsetted_filtered.vcf.gz.tbi "${OUTPATH}/subsetted_filtered.vcf.gz.tbi"
date
EOF

    if command -v ts >/dev/null 2>&1; then
      bash /tmp/RUNNER |& ts
    else
      bash /tmp/RUNNER
    fi

    if [ -n "${DSTAT_PID:-}" ]; then kill "$DSTAT_PID" 2>/dev/null || true; fi
  >>>

  output {
    File usage = "usage.csv"
    File done  = "goodTask"
  }

  runtime {
    docker: docker
    memory: "~{Memory} GB"
    disks: "local-disk ~{Disk} HDD"
    cpu: "~{CPU}"
    preemptible: Preemtible
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
    touch usage.csv
    if command -v dstat >/dev/null 2>&1; then
      dstat -t --cpu --mem --disk --io --freespace > usage.csv 2>&1 & DSTAT_PID=$!
    fi

    cat > /tmp/RUNNER << 'EOF'
set -euo pipefail
eval "$(/bin/micromamba shell hook -s posix)"
micromamba activate base

mkdir -p torun
cd torun

gcloud storage cp -r "~{thisBucket_libfolder}/~{thisLibName}/01_cellSNP/" .
gcloud storage cp    "~{thisBucket_libfolder}/~{thisLibName}/04_VCFSUBSET2/subsetted_filtered.vcf.gz" .
gcloud storage cp    "~{thisBucket_libfolder}/~{thisLibName}/04_VCFSUBSET2/subsetted_filtered.vcf.gz.tbi" .

mkdir -p vireo_final_out
vireo -c 01_cellSNP/ -o vireo_final_out -d subsetted_filtered.vcf.gz -t GT

if [ ! -s "vireo_final_out/summary.tsv" ]; then
  echo "Vireo with doublet detection failed."
  exit 1
else
  echo ok > ../goodTask
fi
ls -l ..

export OUTPATH="~{thisBucket_libfolder}/~{thisLibName}/05_VIREO"
date > .folder
gcloud storage cp .folder "${OUTPATH}/.folder"
gcloud storage cp -r vireo_final_out "${OUTPATH}/"
date
EOF

    if command -v ts >/dev/null 2>&1; then
      bash /tmp/RUNNER |& ts
    else
      bash /tmp/RUNNER
    fi

    if [ -n "${DSTAT_PID:-}" ]; then kill "$DSTAT_PID" 2>/dev/null || true; fi
  >>>

  output {
    File usage = "usage.csv"
    File done  = "goodTask"
  }

  runtime {
    docker: docker
    memory: "~{Memory} GB"
    disks: "local-disk ~{Disk} HDD"
    cpu: "~{CPU}"
    preemptible: Preemtible
    zones: ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
  }
}
