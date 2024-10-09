#!/bin/bash

HOST_IP=$(hostname -i)
PORT_NUM=$(for port in {8700..9000}; do ss -tuln | grep -q ":$port " || { echo "$port"; break; }; done)
echo "***************"
echo "$HOST_IP:$PORT_NUM"
echo "***************"

 trap 'echo signal recieved in BATCH!; kill -15 "${PID}"; wait "${PID}";' SIGINT SIGTERM
trap 'echo signal recieved in BATCH!' SIGINT SIGTERM
podman inspect --format '{{.State.Pid}}' $CONTAINER_ID


IMAGE='us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest'
JUPYTER_CMD="source ~/.profile ; micromamba activate base ; jupyter lab --allow-root --ip='*' --port='8888' --NotebookApp.token='' --NotebookApp.password='' --Session.key='' --no-browser"
CONTAINER_ID=$(podman run --rm --detach -v /broad/macosko/data/discopipeline:/discopipeline -p $PORT_NUM:8888 --entrypoint "/bin/bash" $IMAGE -c "$JUPYTER_CMD")
CONTAINER_PID=$(podman inspect --format '{{.State.Pid}}' $CONTAINER_ID)

wait $CONTAINER_PID
