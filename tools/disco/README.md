### Quick Start guide

1\) Read the background on [slurm](slurm.md)

```
ssh login00.broadinstitute.org
```
```
mkdir -p /broad/macosko/$USER
cd /broad/macosko/$USER
```

2\) Build a podman image:

```
srun --partition=hpcx_macosko --mem 16G --pty bash
podman load -i /broad/macosko/discopipeline/scripts/pipeline-image.tar
wget -O Dockerfile https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/tools/disco/Dockerfile
podman build -f Dockerfile -t myimage --build-arg PASSWORD=$USER .
exit
```

3\) Add these methods to your `.my.bashrc`:

```
IMAGE='myimage'

GET_PORT='HOST_IP=$(hostname -i) ; 
          PORT_NUM=$(for port in {8000..9000}; do ss -an | grep -q :$port || { echo $port; break; }; done) ; 
          echo -e "****************\n$HOST_IP:$PORT_NUM\n****************"'

PODMAN_RUN='exec podman run --rm -it --init --pull never         \
            -v /broad/macosko:/broad/macosko:ro                  \
            -v /broad/macosko/$USER:/broad/macosko/$USER:rw      \
            -p ${PORT_NUM:-}:8787                                \
            -w /broad/macosko/$USER'        

PRINT_TRES='squeue -j $SLURM_JOB_ID -o "%5A %7T %12P %22N %4C %10m %11l"'
DEFAULT_TRES='--partition hpcx_macosko --mem 16G -t 24:00:00'

rstudio() {
    hostname | grep -qv login && echo "Must be on login server" && return 1
    srun --pty -X -C container -J rstudio $DEFAULT_TRES $@ bash -c "$PRINT_TRES ; $GET_PORT ; $PODMAN_RUN $IMAGE rstudio"
}

jupyterlab() {
    hostname | grep -qv login && echo "Must be on login server" && return 1
    srun --pty -X -C container -J jupyterlab $DEFAULT_TRES $@ bash -c "$PRINT_TRES ; $GET_PORT ; $PODMAN_RUN $IMAGE jupyterlab"
}
```

You can start a `tmux` session and run any of these commands:

```rstudio```  

```jupyterlab```  

You can pass in arguments for `srun` after the command, for example:

```rstudio --mem 64G --time 3-00:00:00```

The username is `root` and the password is your username as it appears in `$USER`
