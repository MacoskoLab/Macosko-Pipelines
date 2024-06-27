version 1.0

# Compute disk size and run checks
task getfastqsize {
    input {
        String fastqs_path
        Array[String] samples
        Array[Array[Int]] lanes # [[]] by default (take all lanes of all indexes)
        String memory_multiplier     
        String docker                
    }
    command <<<
        echo "<< starting getdisksize >>"

        # Validate the input
        if ! gsutil ls "~{fastqs_path}" &> /dev/null; then
            echo "ERROR: gsutil ls command failed on input FASTQ path (~{fastqs_path}), does it exist?"
            exit 1
        fi
        if [ ~{length(samples)} -eq 0 ]; then
            echo "ERROR: the list of provided samples is empty, not able to gather FASTQs"
            exit 1
        fi
        if [[ '~{sep="-" samples}' = *" "* ]]; then
            echo "ERROR: sample names are not permitted to have spaces"
            exit 1
        fi

        # Convert from WDL datastructures into bash arrays
        samples=(~{sep=' ' samples})
        lanes=3




        if [ ~{length(lanes)} -eq 1 ] && [ ~{length(lanes[0])} -eq 0 ]; then
            echo "Taking all lanes of all indexes"
        else
            if [ ~{length()} -neq ]
        exit 1



        # Combine the sample name and lanes into a unique id
        if [ ~{length(lanes)} -eq 0 ]; then
            id="~{sample}"
        else
            id="~{sample}_L~{sep='-' lanes}"
        fi
        echo "FASTQ directory: ~{fastq_directory}"
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
            echo "ERROR: size limit reached, increase cap ($(cat SIZE) of 1024 GiB)"
            rm -f SIZE
        fi

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
