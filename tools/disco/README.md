Building an image
-----------------

First, we need to be on a node that supports the container feature. Let's open an interactive session on our partition:

```srun --partition=hpcx_macosko --pty bash```

We've made an example Dockerfile to use for this tutorial. Download this file in your home directory:

```
cd /broad/macosko/$USER
wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/tools/disco/Dockerfile
```

Build the image: (should take about 5 minutes)

```podman build --squash -t myimage .```

Now when we list the built podman images, you should see the newly created image:

```podman images --all```

You'll also see a debian image which was pulled to be used as a base - this can be kept or removed with podman rmi <image>

A few notes on the build before moving on:

`-t myimage`: this is the name to give the image (tag?)
`--squash`: podman stores a copy of the image after each step of the build - this is to make
If you rebuild, the new image does not overwrite the old one. get the "latest" tag and the previous image become untagged. They take up space so Make sure to stop all containers using the old image and podman rmi it

Podman images are stored in `/local/podman/containers-user-$UID`, and this folder is NOT shared across the various nodes. This means that you won't be able to make containers using this image on a different node. We can run a command (TODO)

Exit the srun session - this step is complete

```exit```

These commands do not have to be run again, unless you decide to update the Dockerfile. The image is now available anytime to be used as a container template.

Intro to containers
-------------------

Step 1: read the background on slurm `here`
Step 2: create a podman image `here`