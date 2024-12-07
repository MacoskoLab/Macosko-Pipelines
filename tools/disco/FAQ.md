Frequently Asked Questions

---

Q: `srun: command not found`

A: Make sure you are on login00.broadinstitute.org

---

Q: When I run `srun` it hangs forever or gives this error:

```
srun: error: Unable to create step for job: Requested node configuration is not available
```

A: `srun` will hang forever if you are already in a node - make sure you are on login00.broadinstitute.org. Alternatively, the node may be busy. Try requesting requesting fewer resources or monitoring usage with TODO

---

Q: `podman: command not found`

A: Make sure you are on a node with container support - this can be specified with `-C container`

---

Q: When I call `podman run` on "myimage", this happens:

```
Error: myimage: image not known
```
or
```
? Please select an image: 
  ▸ registry.access.redhat.com/myimage:latest
    registry.redhat.io/myimage:latest
    docker.io/library/myimage:latest
```
or
```
Resolving "myimage" using unqualified-search registries (/etc/containers/registries.conf)
Trying to pull registry.access.redhat.com/myimage:latest...
Trying to pull registry.redhat.io/myimage:latest...
Trying to pull docker.io/library/myimage:latest...
Error: 3 errors occurred while pulling:
 * initializing source docker://registry.access.redhat.com/myimage:latest: reading manifest latest in registry.access.redhat.com/myimage: name unknown: Repo not found
 * initializing source docker://registry.redhat.io/myimage:latest: unable to retrieve auth token: invalid username/password: unauthorized: Please login to the Red Hat Registry using your Customer Portal credentials. Further instructions can be found here: https://access.redhat.com/RegistryAuthentication
 * initializing source docker://myimage:latest: reading manifest latest in docker.io/library/myimage: requested access to the resource is denied
```

A: Podman could not find "myimage". Run `podman images -a` to verify that "myimage" exists. Podman images are not shared between nodes, so try making a copy or going to the node where the image was created.

---

Q: How can I log podman into GCP?

1. `ssh login00.broadinstitute.org`
2. `use Google-Cloud-SDK`
3. srun/ssh into the podman node
4. `gcloud auth print-access-token --quiet | podman login -u oauth2accesstoken --password-stdin us-central1-docker.pkg.dev`

---

Q: What are some helpful commands for managing podman?

A:
* `podman login`: Log in to a container registry (e.g. our GCP artifact registry)
* `podman system df`: Show podman disk usage
* `podman system prune`: Remove all unused pods, containers, images, networks, and volume data

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

Alternatively, you could just SCP it directly:

```
srun --partition=hpcx_macosko --mem 20G --pty bash
export TMPDIR="/broad/macosko/$USER"

echo 'Adding connections:'
for NODE in $(sinfo -r -h -o "%n %f" | awk '$2 ~ /container/ {print $1}') ; do
    echo $NODE
    podman system connection add --identity /home/unix/$USER/.ssh/id_rsa $NODE $USER@$NODE
done
podman system connection default $(hostname -s)

IMAGE='myimage'

echo 'Distributing the image:'
for NODE in $(sinfo -r -h -o "%n" | grep slurm-bits-d) ; do
    echo $NODE
    podman image scp "$(hostname -s)::$IMAGE" "$NODE::"
done

exit
```

---

Q: What is the easiest way to transfer files in/out of a container?

Option 1: use `-v` to share a folder inside the container with a folder outside the container

Option 2: use `podman cp` to copy files/folders between a container and the local filesystem

Option 3: tar?

Option 4: mount the container's entire filesystem in a user namespace shell:

```
podman unshare
podman mount mycontainer
```

You should see a directory printed to the terminal which contains the root file system of the container - you can now copy in/out any data that you need. Once finished, exit the shell:

```exit```

---

Q: Can I back up a container?

Option 1: `podman commit` creates an image based on a changed container

Option 2: `podman export` exports the filesystem of a container

Option 3: `podman container clone` creates a copy of a container, recreating the original with an identical configuration

---

Q: Can I change port bindings or volume mounts after creating the container?

A: No. You need to create a new container 

---

Help! The `srun` session ended while a container was still running and now the podman state is corrupted. How to fix?

Step 1: Log out of the node and log back in 

Step 2: Run `podman system migrate`

Step 3: Run `podman ps --all --external` - you should see the container has the status "Stopping"

Step 4: The container cannot be repaired - we will instead copy its state to a new image and use that image to generate a new container. The contents of memory will be lost, but everything on disk will be preserved. Run the command `podman commit mycontainer newimage`

Step 4: The container cannot be repaired - to delete it, first make sure that there are no other containers running on the node, then run these commands:

```
rm -rf /tmp/containers-user-$UID
rm -rf /tmp/podman-run-$UID
```

You should now be able to `podman rm -f <NAME>` the orphaned container

Note: When the time limit of `srun` elapses, the container will be sent a SIGTERM signal which should cleanly clean up its state. These instructions will instead usually have to be run anytime the terminal session containing the `srun` is prematurely closed.

Step 5: You can now use "newimage" to make a container

To remove the broken container, you must remove the image the container was made from using `podman rmi --force myimage`. Make sure all containers from this image have been backed up as they will be deleted too.

---

Q: Can I completely reset podman state?

A: Yes. Run these commands on each node that has podman state:

```
rm -rf /tmp/containers-user-$UID
rm -rf /tmp/podman-run-$UID
podman rmi -a -f
echo y | podman system reset
```

Any containers found left in `podman ps --all --external` will have to be manually removed, then run `podman system reset` again
