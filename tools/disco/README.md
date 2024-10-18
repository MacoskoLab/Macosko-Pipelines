Introduction to login00
=======================

The login server `login00.broadinstitute.org` has SLURM enabled. We can use it to submit jobs to various nodes
 nodes that jobs can be submitted to
 
You can login to this 
``` ssh ```

Introduction to nodes
=====================

There are a total of three partitions, each with their own set of nodes
* disco
* gpu
* hpcx_macosko

The disco partition has 12 nodes:
* slurm-bits-d001
* slurm-bits-d002
* slurm-bits-d003
* slurm-bits-d004
* slurm-bits-d005
* slurm-bits-d006
* slurm-bits-rh7-d001
* slurm-bits-rh7-d002
* slurm-bits-rh7-d003
* slurm-bits-rh7-d004
* slurm-bits-rh7-d005
* slurm-bits-rh7-d006
Available Features:
`slurm-bits-d[001-006]`: RedHat8, container
`slurm-bits-rh7-d[001-006]`: RedHat7, dotkit
Each node has 128 cores and 755GiB RAM

The gpu partition has 3 GPU nodes:
* slurm-gpu-d001
* slurm-gpu-d002
* slurm-gpu-d003
Available Features:
`slurm-gpu-d[001-003]`: RedHat8, container
Each GPU node has 64 cores, 502GiB RAM, and 4 NVIDIA A30 with 24GiB VRAM

Finally, our lab has exclusive access to one giant node on the hpcx_macosko partition:
`slurm-bits-bigmem-d002`: RedHat8, container
This node has 64 cores and 4030GiB of RAM

Type `sinfo` or `scontrol show nodes` to see this information in more detail

`container` means the node has podman, apptainer, and singularity
`dotkit` means the node has the familiar list of "use" modules (e.g. Google-Cloud-SDK)

Introduction to SLURM
=====================

There are two ways to run jobs:

`srun` runs the job interactively in the terminal. Output is printed in real-time and the job can be control-C'd to end. Closing the terminal will also end the job.
`sbatch` submits the job to be handled by slurm. It runs in the background and writes output to a log file. The terminal can be safely closed. Analogous to qsub

The rest of this tutorial will use `srun` due to its interactive

This command will put you in a bash terminal on an available node:
`srun --pty bash`

TODO: explain what this does, and the default parameters, mounted files

By default, these are the resources:
* --partition = disco
* --time = 01:00:00 (1 hour)
* -c --cpus-per-task = 1
* --mem = 1G per CPU

Here are some helpful commands:

`squeue -u $USER` lists all your jobs
`squeue -t RUNNING` lists all running jobs
`scontrol show job <ID>` lists information about a job
`scancel <ID>` cancels a job
`scancel -u $USER` cancels all your jobs

Any modules loaded in the login server will propagate to the srun node

Introduction to podman
======================

There is no way to get root access on the Broad HPC cluster. Thus, all commands that required `sudo` such as `apt install` are unavailable. As such, we will instead be setting up our environment inside of a rootless podman container.

Here are three important terms:
* Dockerfile: a recipe to build a single image
* Image: a read-only template used for creating containers
* Container: an isolated environment in which we can install packages and run code

We have a sample Dockerfile available. In the next section we'll show how to build an image using this recipe, and in the section after we'll use this image to create runnable containers.

Building an image
=================

The first thing we need to do is make a podman image. To do this, we need to be on a node that supports the container feature. Let's open an interactive session on our partition:

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

Podman images are stored in `/local/podman/containers-user-$UID`, and this folder is NOT shared across the various nodes. We can run a command (TODO)

Exit the srun session - this step is complete

```exit```

These commands do not have to be run again, unless you decide to update the Dockerfile. The image is now available anytime to be used as a container template.

Intro to containers
===================
