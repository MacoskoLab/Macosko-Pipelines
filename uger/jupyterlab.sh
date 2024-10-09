#!/bin/bash

HOST_IP=$(hostname -i)
PORT_NUM=$(for port in {8700..9000}; do ss -tuln | grep -q ":$port " || { echo "$port"; break; }; done)
echo "***************"
echo "$HOST_IP:$PORT_NUM"
echo "***************"

trap "echo TRAPTRIGGERED ; podman stop jupyterlab$PORT_NUM ; exit" SIGTERM SIGINT

IMAGE='us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest'
JUPYTER_CMD="source ~/.profile ; micromamba activate base ; jupyter lab -y --allow-root --ip='*' --port='8888' --NotebookApp.token='' --NotebookApp.password='' --Session.key='' --no-browser"
podman run --rm --name jupyterlab$PORT_NUM -v /broad/macosko/data/discopipeline:/discopipeline -p $PORT_NUM:8888 --entrypoint "/bin/bash" $IMAGE -c "$JUPYTER_CMD"
