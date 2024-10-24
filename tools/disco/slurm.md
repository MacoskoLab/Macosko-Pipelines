Introduction to login00
-----------------------

The login server `login00.broadinstitute.org` has SLURM enabled. We can use it to submit jobs to various nodes that have large amounts of computing resources.
 
You can login to this server using ssh:

``` ssh login00.broadinstitute.org ```

To make logging in easier, run these commands on your local machine:

```
ssh-keygen -t rsa -b 4096 -f ~/.ssh/uger -N ""
ssh-copy-id -i ~/.ssh/uger.pub $USER@login00.broadinstitute.org
```
```
echo -e "\nHost l0\n\tHostName login00.broadinstitute.org\n\tUser $USER\n\tIdentityFile ~/.ssh/uger" >> ~/.ssh/config
```

Now you should be able run `ssh l0` to log in without needing a password

Introduction to nodes
---------------------

There are a total of three partitions, each with their own set of nodes:
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
---------------------

There are two ways to run jobs:

`srun` runs the job interactively in the terminal. Output is printed in real-time and the job can be control-C'd to end. Closing the terminal will also end the job.
`sbatch` submits the job to be handled by slurm. It runs in the background and writes output to a log file. The terminal can be safely closed. Analogous to qsub

The rest of this tutorial will use `srun` due to its interactive nature

This command will put you in a bash terminal on an available node:

```srun --pty bash```

By default, these are the arguments:
* --partition = disco
* -t --time = 01:00:00 (1 hour)
* -c --cpus-per-task = 1
* --mem = 1G (per CPU)

Here are some helpful commands:
* `squeue -u $USER`: lists all your jobs
* `squeue -t RUNNING`: lists all running jobs
* `scontrol show job <ID>`: lists information about a job
* `scancel <ID>`: cancels a job
* `scancel -u $USER`: cancels all your jobs

Any modules loaded in the login server will propagate to the srun node

Introduction to podman
----------------------

Commands requiring root permisison cannot be run on the Broad HPC cluster. This means that all commands which require `sudo` such as `apt install` are unavailable. Therefore, we will instead be setting up our environment inside of a podman container.

Inside a container, we have root access. We can install packages, configure an environment, and run code.

A container is sealed off from the surrounding environment.
Inside has its own filesystem, process tree, and network stack. 

Here are three important terms:
* Dockerfile: a recipe to build a single image
* Image: a read-only template used for creating containers
* Container: an isolated environment in which we can install packages and run code

In the next section we'll show how to build an image using a sample Dockerfile, and in the section after we'll use this image to create runnable containers.

Helpful links
-------------
https://broad.service-now.com/kb_view.do?sys_kb_id=a6c74cb147d6a51411484438946d430e
https://backstage.broadinstitute.org/docs/default/component/disco-docs
