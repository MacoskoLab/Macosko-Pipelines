## Dotkit
`use -l`: List all available dotkits  
`use -l <text>`: Search available dotkits

Example dotkits:
* `use UGER`
* `use Anaconda3`
* `use Google-Cloud-SDK`


## Helpful commands
`qhost`: Show the status of UGER hosts  
`qstat`: Show the status of UGER jobs  
`qstat -j <jobID>`: Print job information  
`qstat -u “*”`: Compare the priority of your job to others in UGER  
`qacct -j <jobID>`: List completed job details  
`qdel <jobID>`: Delete an unwanted job  
`qdel -f <jobID>`: Force delete stuck (dr and dt) jobs  


## Interactive jobs
`/broad/software/scripts/ish` is a wrapper for `qrsh -q interactive -now n -l h_rt=36:00:00 $@`

Example:  
```/broad/software/scripts/ish -l h_vmem=8g -pe smp 3 -binding linear:3 -l h_rt=18:00:0 -P macosko_lab -l os='RedHat7'```


## Batch jobs

Option 1: `qsub -b y <binary>`

Option 2: `qsub <script.sh>`

`#$` lines indicate flags for qsub

Example:
```
#! /bin/bash
#$ -N FirstScript

source /broad/software/scripts/useuse
reuse -q BLAST

blastall <params>
```

Command line flags override job script settings


## Job submission defaults
|                     |                                               |
|---------------------|-----------------------------------------------|
| queue               | “short” (has a two hour runtime limit)        |
| memory              | 1G                                            |
| cores               | 1
| execution directory | your home directory                           |
| output filename     | `<jobName>.o<jobID>` (in execution directory) |
| stderr filename     | `<jobName>.e<jobID>` (in execution directory) |
| name                | “word” in job command                         |


## Job submission options
`-N name`: The name of the job  
`-o path`: Standard output  
`-e path`: Standard error  
`-j y|n`: Merge standard error into standard output  
`-cwd` or `-wd working_dir`: Job exec location  
`-w e`: Jobs with invalid requests will be rejected  
`-pe smp n -binding linear:n`: Bind the job on `n` successive cores  
`-l resource=value`: Launch the job in a UGER queue meeting the given resource request list
* `-l h_rt=<hh:mm:ss>`: Hard run time limit
* `-l h_vmem=<size>`: Hard virtual memory limit (per processor slot)


## UGER host resources:
`os`: RedHat7 (49), RedHat8 (21)  
`gpu`: 1 (3)  
`processor`: Intel(R) Xeon(R) Platinum 8358 CPU @ 2.60GHz (70)  
`arch`: lx-amd64 (70)  
`num_proc`: 64 (30), 128 (40)  
`mem_total`: 755G (26), 770G (1), 251G (21), 502G (3), 250G (7), 754G (11), 4030G (1)  


## Helpful links:
[Intro to UGER - Nov 2016](https://data.broadinstitute.org/bits_demo/user_education_sessions/Intro2UGER/I2U-Nov2106/Intro2UGER.pdf)
[How to submit a job using qsub](https://bioinformatics.mdc-berlin.de/intro2UnixandSGE/sun_grid_engine_for_beginners/how_to_submit_a_job_using_qsub.html)
