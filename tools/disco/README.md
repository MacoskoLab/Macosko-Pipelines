1\) Read the background on [slurm](slurm.md)

2\) Build a podman [image](images.md)

3\) Add these methods to your `.my.bashrc`:

```
IMAGE='myimage'

GET_PORT='HOST_IP=$(hostname -i) ; 
          PORT_NUM=$(for port in {8000..9000}; do ss -an | grep -q :$port || { echo $port; break; }; done) ; 
          echo -e "****************\n$HOST_IP:$PORT_NUM\n****************"'

PODMAN_RUN='exec podman run --rm -it --init                      \
            -v /broad/macosko:/broad/macosko:ro                  \
            -v /broad/macosko_storage:/broad/macosko_storage:ro  \
            -v /broad/macosko/$USER:/broad/macosko/$USER:rw      \
            -p ${PORT_NUM:-}:8787                                \
            -w /broad/macosko/$USER'        

rstudio() {
    hostname | grep -qv login && echo "Must be on login server" && exit 1
    srun --pty -X -C container -J rstudio -t 24:00:00 $@ bash -c "$GET_PORT ; $PODMAN_RUN $IMAGE rstudio"
}

jupyterlab() {
    hostname | grep -qv login && echo "Must be on login server" && exit 1
    srun --pty -X -C container -J jupyterlab -t 24:00:00 $@ bash -c "$GET_PORT ; $PODMAN_RUN $IMAGE jupyterlab"
}

```

You can start a `tmux` session and run any of these commands:

```rstudio```

```jupyterlab```

The username is `root` and the password is your username as it appears in `$USER`
