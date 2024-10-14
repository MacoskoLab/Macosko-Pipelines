#!/bin/bash

HOST_IP=$(hostname -i)
PORT_NUM=$(for port in {8700..9000}; do ss -an | grep -q ":$port " || { echo "$port"; break; }; done)
echo "***************"
echo "$HOST_IP:$PORT_NUM"
echo "***************"

IMAGE='us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest'
JUPYTER_CMD="source ~/.profile ; micromamba activate base ; jupyter lab --allow-root --ip='*' --port='8888' --NotebookApp.token='' --NotebookApp.password='' --Session.key='' --no-browser"
exec podman run --rm --name "jupyterlab$PORT_NUM" -v /broad/macosko/data/discopipeline:/discopipeline -p $PORT_NUM:8888 --entrypoint "/bin/bash" $IMAGE -c "$JUPYTER_CMD"

# Add this function to your .my.bashrc
jupyterlab() {
    srun -X -C container -J jupyterlab -t 30-00:00:00 "$@" /broad/macosko/data/discopipeline/scripts/jupyterlab.sh
}
