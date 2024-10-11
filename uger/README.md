### Podman commands

Logging in on all avaiable nodes (run in login00.broadinstitute.org)
```
use Google-Cloud-SDK
LOGIN_CMD="gcloud auth print-access-token --quiet | podman login -u oauth2accesstoken --password-stdin us-central1-docker.pkg.dev"

PARTITION_NODE_LIST=$(sinfo -N -h -o "%P %N" | sed 's/\*//')
while read -r partition node <&3; do
  echo "Logging in on $partition $node:"
  srun --partition=$partition --nodelist=$node bash -c "$LOGIN_CMD"
done 3<<< "$PARTITION_NODE_LIST"
```

