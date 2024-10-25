Building an image
-----------------

First, we need to be on a node that supports the container feature. Let's open an interactive session on our partition:

```srun --partition=hpcx_macosko --mem=10G --pty bash```

We've made an example Dockerfile to use for this tutorial. Download this file in your home directory:

```
cd /broad/macosko/$USER
wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/tools/disco/Dockerfile
```

Build an image using this recipe: (should take about 5 minutes)

```podman build --squash -t myimage .```

List the built podman images:

```podman images -a```

You should see a newly created image called "myimage"

Podman images are stored in `/local/podman/containers-user-$UID`, and this folder is NOT shared across the various nodes. This means that you currently won't be able to make containers using this image on a different node. See the FAQ for a discussion on how to send a copy of the image to other nodes.

Exit the srun session - this step is complete

```exit```

These commands do not have to be run again, unless you decide to update the Dockerfile. The image is now available on this node anytime to be used as a container template.

Summary
-------
Log into `login00.broadinstitute.org` and run this command:
```
srun --partition=hpcx_macosko --pty bash
wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/tools/disco/Dockerfile
podman build -t myimage  --build-arg PASSWORD=$USER .
exit
```
