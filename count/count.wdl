version 1.0

workflow count {
    String pipeline_version = "1.0.0"
    input {
        String bcl
        String samplesheet
        String count_output = "gs://"+bucket+"/counts/"+basename(bcl,"/")
        String log_output = "gs://"+bucket+"/logs/"+basename(bcl,"/")
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-bcl2fastq/img:latest"
    }
    parameter_meta {
    }
    
    output {
    }
}
