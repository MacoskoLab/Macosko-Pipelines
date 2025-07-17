version 1.0

task mkfastq {
    input {
        String bcl
        String samplesheet
        Int mem_GB
        Int disk_GB
        String params
        String docker
    }
    command <<<
    set -euo pipefail

    BUCKET="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
    fastq_dir="gs://$BUCKET/fastqs/$(basename ~{bcl})"

    echo "==================== START BCL-CONVERT ===================="

    # Assert output does not already exist in bucket
    if gsutil -q stat "$fastq_dir/Logs/FastqComplete.txt" ; then
        echo "Error: Output already exists in the bucket!"
        exit 1
    fi

    echo "----- Downloading BCL -----"
    gcloud storage cp -r "~{bcl}" .
    mv $(basename ~{bcl}) BCL
    ls -1 BCL

    echo "----- Downloading SampleSheet -----"
    gcloud storage cp "~{samplesheet}" SampleSheet.csv
    cat SampleSheet.csv
    
    echo "----- Running bcl-convert -----"
    bcl-convert                        \
        --bcl-input-directory=BCL      \
        --output-directory=fastqs      \
        --sample-sheet=SampleSheet.csv \
        ~{params}

    rm fastqs/Undetermined_S0_*

    echo "----- Uploading results -----"
    gcloud storage cp -r fastqs/* "$fastq_dir/"

    # Assert output had been successfully uploaded to the bucket
    if ! gsutil -q stat "$fastq_dir/Logs/FastqComplete.txt" ; then
        echo "Error: Output not found in the bucket!"
        exit 1
    fi

    echo "==================== END BCL-CONVERT ===================="

    >>>
    runtime {
        cpu: 8
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 0
    }
}


task getdisksize {
    input {
        String bcl
        String docker
    }
    command <<<
        echo "----- START CHECK -----"

        # Assert that the BCL exists
        if ! gsutil ls "~{bcl}" &> /dev/null; then
            echo "ERROR: gsutil ls failed on ~{bcl}"
            exit 1
        fi

        # Get the size of the BCL * 2.5
        gsutil du -sc "~{bcl}" | grep total | 
        awk '{size=$1/1000/1000/1000 ; size=size*2.5 ; if (size<128) size=128 ; printf "%d\n", size}' > SIZE

        # Assert that the disksize file exists
        if [ ! -s SIZE ]; then
            echo "ERROR: gsutil du failed on ~{bcl}"
            exit 1
        fi

        echo "----- END CHECK -----"
    >>>
    output {
        Int disk_GB = read_int("SIZE")
    }
    runtime {
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 8 HDD"
        docker: docker
        preemptible: 0
  }
}


workflow bcl_convert {
    input {
        String bcl
        String samplesheet
        String params = "--strict-mode=true"
        String docker = "us-central1-docker.pkg.dev/velina-208320/terra/software-image:latest"
    }
    call getdisksize {
        input:
            bcl = bcl,
            docker = docker
    }
    call mkfastq {
        input:
            bcl = bcl,
            samplesheet = samplesheet,
            mem_GB = 64,
            disk_GB = getdisksize.disk_GB,
            params = params,
            docker = docker
    }
}