Frequently Asked Questions

---

Q: `srun: command not found`

A: Make sure you are on login00.broadinstitute.org

---

Q: When I run `srun` it hangs forever - why?

A: `srun` will hang forever if you are already in a node - make sure you are on login00.broadinstitute.org. Alternatively, the node may be busy. Try requesting requesting fewer resources or monitoring usage with TODO

---

Q: When I call `podman run` on "myimage", this happens:

```
? Please select an image: 
  ▸ registry.access.redhat.com/myimage:latest
    registry.redhat.io/myimage:latest
    docker.io/library/myimage:latest
```

A: Podman could not find "myimage". Run `podman images -a` to verify that "myimage" exists. Podman images are not shared between nodes, so try making a copy or going to the node where the image was created.

---

Q: Can I share images between nodes?

A: Yes. First same the image to a .tar file:

```podman save -o myimage.tar myimage```

Then run this commmand to distribute it across all nodes:

```
COMMAND="podman load -i myimage.tar"
COMMAND="if command -v podman &> /dev/null; then $COMMAND; else echo podman not found; fi"
PARTITION_NODE_LIST=$(sinfo -N -h -o "%P %N" | sed 's/\*//')
while read -r partition node <&3; do
  echo "$partition $node:"
  srun --partition=$partition --nodelist=$node bash -c "$COMMAND"
done 3<<< "$PARTITION_NODE_LIST"
```

Alternatively, you could just run this command to scp it directly:

```
podman image scp
```

---

Help! The `srun` session ended while "mycontainer" was still running and now the podman state is corrupted. How to fix?

Step 1: Log out of the node and log back in 

Step 2: Run `podman system migrate`

Step 3: Run `podman ps --all --external` - you should see the container has the status "Stopping"

Step 4: The container cannot be repaired - we will instead copy its state to a new image and use that image to generate a new container. The contents of memory will be lost, but everything on disk will be preserved. Run the command `podman commit mycontainer newimage`

Step 5: You can now use "newimage" to make a container

To remove the broken container, you must remove the image the container was made from using `podman rmi --force myimage`. Make sure all containers from this image have been backed up as they will be deleted too.

---

Q: Can I completely reset podman state?

A: Yes. Run these commands on each node that has podman state:

```
rm -rf /tmp/containers-user-$UID
rm -rf /tmp/podman-run-$UID
podman system reset
```
