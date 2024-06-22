# Compute disk size needed for cellranger count and run checks
task getfastqsize {
    input {
        String fastqs
        String sample
        Array[Int] lanes
        String output_path
        String log_output_path       
        String docker                
    }
    command <<<
        echo "<< starting getdisksize >>"

        # Combine the sample name and lanes into a unique id
        if [ ~{length(lanes)} -eq 0 ]; then
            id="~{sample}"
        else
            id="~{sample}_L~{sep='-' lanes}"
        fi
        echo "FASTQ path: ~{fastqs}"
        echo "sample: ~{sample}"
        echo "lanes: ~{sep=',' lanes}"
        echo "id: $id" ; echo
        echo $id > ID

        # Get the fastq files and their total size
        gsutil ls -r ~{fastqs} | fgrep ".fastq.gz" | fgrep "~{sample}_S" | fgrep -v "_I1_" | fgrep -v "_I2_" > PATHS
        if [ ~{length(lanes)} -gt 0 ]; then
            lanes=(~{sep=' ' lanes})
            regex=$(printf "_L0*%s_ " "${lanes[@]}" | sed 's/ $//' | sed 's/ /|/g')
            grep -P $regex PATHS > temp_file && mv temp_file PATHS
        fi
        cat PATHS | xargs gsutil du -sc | grep total | awk '{size=$1/1024/1024/1024 ; size=size*6+20 ; if (size<128) size=128 ; printf "%d\n", size}' > SIZE
        if [ ~{length(lanes)} -gt 0 ]; then
            lanes=(~{sep=' ' lanes})
            for lane in "${lanes[@]}"; do
                if ! grep -q "_L0*${lane}_" PATHS; then
                    echo "ERROR: No fastqs found for lane ${lane}"
                    rm -f SIZE
                fi
            done
        fi

        echo "FASTQ files:"
        cat PATHS; echo

        echo "Total size (GiB):"
        cat SIZE; echo

        # Assert that the fastqs exist
        if [[ ! -s PATHS ]]
        then
            echo "ERROR: gsutil ls command returned a blank file"
            rm -f SIZE
        fi
        if [[ ! -s SIZE ]]
        then
            echo "ERROR: gsutil du command failed on input fastqs"
            rm -f SIZE
        fi

        # Assert that the disksize is not too large
        if [[ $(cat SIZE) -gt 6000 ]]
        then
            echo "ERROR: cellranger-count disk size limit reached, increase cap"
            rm -f SIZE
        fi

        # Assert that the count output is blank (avoid overwiting)
        output_path="~{output_path}"
        if gsutil ls "${output_path%/}/$id" &> /dev/null
        then
            echo "ERROR: cellranger-count output already exists"
            rm -f SIZE
        else
            echo "Output path: ${output_path%/}/$id"
        fi

        # Assert that the paths are actually gs:// paths
        [[ ! "~{fastqs}" =~ gs:// ]] && echo "ERROR: fastq_path does not contain gs://" && rm SIZE
        [[ ! "~{output_path}" =~ gs:// ]] && echo "ERROR: output_path does not contain gs://" && rm SIZE
        [[ ! "~{log_output_path}"   =~ gs:// ]] && echo "ERROR: log_output_path does not contain gs://" && rm SIZE

        echo "<< completed getdisksize >>"
    >>>
    output {
        String id = read_string("ID")
        Int disksize = read_int("SIZE")
        Array[String] fastq_paths = read_lines("PATHS")
    }
    runtime {
        docker: docker
        memory: "4 GB"
        disks: "local-disk 8 HDD"
        cpu: 1
        preemptible: 0
  }
}
