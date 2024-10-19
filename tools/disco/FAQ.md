```srun: command not found```

make sure you are on login00.broadinstitute.org

```
podman run --rm -it myimage bash
? Please select an image: 
  ▸ registry.access.redhat.com/myimage:latest
    registry.redhat.io/myimage:latest
    docker.io/library/myimage:latest
```
podman could not find the image "myimage" when making the container. Run `podman images -a` to verify that "myimage" exists. The podman images are not shared between node, so if you don't see it there try going to the node where the image was created.
