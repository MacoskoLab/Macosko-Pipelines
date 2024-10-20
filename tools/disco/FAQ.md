```srun: command not found```

Make sure you are on login00.broadinstitute.org

---

When I run `srun` it hangs forever - why?

If you are already on a node `srun` will hang forever - make sure you are on login00.broadinstitute.org.

Alternatively, the node may be busy. Try requesting requesting fewer resources or monitoring usage with TODO

---

When I call `podman run` on "myimage", this happens:

```
? Please select an image: 
  ▸ registry.access.redhat.com/myimage:latest
    registry.redhat.io/myimage:latest
    docker.io/library/myimage:latest
```
podman could not find the image "myimage". Run `podman images -a` to verify that "myimage" exists. podman images are not shared between nodes, so try going to the node where the image was created.

---

Can I move images between nodes?

yes

---

Help! The `srun` session ended while a container was still running and now the podman state is corrupted. How to fix?

Step 1: Log out of the node and log back in 

Step 2: Run `podman system migrate`

Step 3: Run `podman ps --all --external` - you should see the container has the status "Stopping"

Step 4: The container cannot be repaired - we will instead copy its state to a new image and use that image to generate a new container. The contents of memory will be lost, but everything on disk will be preserved. Run the command `podman commit mycontainer newimage`

Step 5: You can now use this image to make a new container
