```srun: command not found```

make sure you are on login00.broadinstitute.org

---

```
podman run --rm -it myimage bash
? Please select an image: 
  ▸ registry.access.redhat.com/myimage:latest
    registry.redhat.io/myimage:latest
    docker.io/library/myimage:latest
```
podman could not find the image "myimage". Run `podman images -a` to verify that "myimage" exists. podman images are not shared between nodes, so if you don't see it try going to the node where the image was created.

---

Can I move images between nodes?

yes
