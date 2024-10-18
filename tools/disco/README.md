# Introduction to Disco

The login server `login00.broadinstitute.org` has SLURM enabled.

Currently there are 12 different nodes on the disco partition that jobs can be submitted to:
```
slurm-bits-d001
slurm-bits-d002
slurm-bits-d003
slurm-bits-d004
slurm-bits-d005
slurm-bits-d006
slurm-bits-rh7-d001
slurm-bits-rh7-d002
slurm-bits-rh7-d003
slurm-bits-rh7-d004
slurm-bits-rh7-d005
slurm-bits-rh7-d006
```
Available Features:
`slurm-bits-d[001-006]`: RedHat8, container
`slurm-bits-rh7-d[001-006]`: RedHat7, dotkit
Each node has 128 cores and 755GiB RAM

There are also a gpu partition with 3 GPU nodes:
```
slurm-gpu-d001
slurm-gpu-d002
slurm-gpu-d003
```
Available Features:
`slurm-gpu-d[001-003]`: RedHat8, container
Each GPU node has 64 cores, 502GiB RAM, and 4 NVIDIA A30 with 24GiB VRAM

Finally, our lab has exclusive access to one giant node on the hpcx_macosko partition:
`slurm-bits-bigmem-d002`: RedHat8, container
This node has 64 cores and 4030GiB of RAM

Type `sinfo` or `scontrol show nodes` to see this information in more detail

`container` means the node has podman, apptainer, and singularity
`dotkit` means the node has the familiar list of "use" modules (e.g. Google-Cloud-SDK)
