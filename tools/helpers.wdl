version 1.0

# Compute disk size and run checks
task getfastqsize {
    input {
        String fastqs
        String sample
        Array[Int] lanes
        String memory_multiplier
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
        cat PATHS | xargs gsutil du -sc | grep total | awk '{size=$1/1024/1024/1024 ; size=size*~{memory_multiplier} ; if (size<96) size=96 ; printf "%d\n", size}' > SIZE
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
        if [[ $(cat SIZE) -gt 1024 ]]
        then
            echo "ERROR: size limit reached, increase cap (1024 GiB)"
            rm -f SIZE
        fi

        # Assert that the paths are actually gs:// paths
        [[ ! "~{fastqs}" =~ gs:// ]] && echo "ERROR: fastq_path does not contain gs://" && rm -f SIZE
        [[ ! "~{output_path}" =~ gs:// ]] && echo "ERROR: output_path does not contain gs://" && rm -f SIZE
        [[ ! "~{log_output_path}"   =~ gs:// ]] && echo "ERROR: log_output_path does not contain gs://" && rm -f SIZE

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
