## Podman commands
Log in:
```
use Google-Cloud-SDK
COMMAND="gcloud auth print-access-token --quiet | podman login -u oauth2accesstoken --password-stdin us-central1-docker.pkg.dev"
```
List the images:
```
COMMAND="podman images --all"
```
List the containers:
```
COMMAND="podman ps --all --external"
```
Run the command on all nodes (run in login00.broadinstitute.org):
```
COMMAND="if command -v podman &> /dev/null; then $COMMAND; else echo podman not found; fi"
PARTITION_NODE_LIST=$(sinfo -N -h -o "%P %N" | sed 's/\*//')
while read -r partition node <&3; do
  echo "$partition $node:"
  srun --partition=$partition --nodelist=$node bash -c "$COMMAND"
done 3<<< "$PARTITION_NODE_LIST"
```
Helpful commands:
Pull image: ```podman pull <IMAGE>```

Delete container: ```podman rm <ID>```

Delete image: ```podman rmi <ID>```

```podman system migrate```

```podman system reset```

```rm -rf /tmp/containers-user-$UID``` (WARNING)

```rm -rf /tmp/podman-run-$UID``` (WARNING)

## Docker images
### us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest
R, Python, Julia, JupyterLab
