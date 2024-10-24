1\) Read the background on [slurm](slurm.md)

2\) Build a podman [image](images.md)

```
export IMAGE='myimage'
export JUPYTER_CMD='source ~/.profile ; micromamba activate base ; jupyter lab --allow-root --ip="*" --port="8888" --NotebookApp.token="" --NotebookApp.password="" --Session.key="" --no-browser'
PODMAN_CMD='HOST_IP=$(hostname -i) ; 
            PORT_NUM=$(for port in {8000..9000}; do ss -an | grep -q ":$port " || { echo "$port"; break; }; done) ; 
            echo "***************" ; 
            echo "$HOST_IP:$PORT_NUM" ; 
            echo "***************" ; 
            exec podman run --rm --init --name jupyterlab-$HOST_IP:$PORT_NUM 
            -v /broad/macosko:/broad/macosko:ro 
            -v /broad/macosko_storage:/broad/macosko_storage:ro 
            -v /broad/macosko/$USER:/broad/macosko/$USER:rw 
            -p $PORT_NUM:8888 
            $IMAGE bash -c "$JUPYTER_CMD"'

srun -X -C container -J jupyterlab -t 10:00:00 $@ bash -c "$PODMAN_CMD"
```
