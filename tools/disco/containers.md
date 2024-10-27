Intro to containers
-------------------

Let's start a new srun session, this time giving us more time and memory:

srun --partition=hpcx_macosko --time 1-00:00:00 --mem 32G --pty bash

Finally, it's time to make a container. This command creates a new container:

podman create --name mycontainer --detach --init myimage

--init enables the container to be run in the background

You can list running podman containers with this command:

podman ps -a

You should see a new container created. Now let's attach to it:

podman exec -it mycontainer /bin/bash

You are now running a bash terminal inside the container! By default, you will be root with home directory /root. If you take a look around with ls you'll see that there is a whole new filesystem in here. Exit anytime:

exit

If you podman ps you'll see that the container is still running. Run podman stop to halt the execution:

podman stop -t 0 mycontainer

Stopped containers can be seen with podman ps --all. You can use podman start to restart the stopped container:

podman start mycontainer

Any changes made to the container are persistent and will remain after restarting the container. . however, just like images, only this node has them. containers are stored in /local and /tmp/containers-user-$UID. Running containers /tmp/podman-run-$UID (running containers)

podman stop -t 0 mycontainer
podman rm mycontainer

It is critical to ensure the container is stopped before the srun session ends. Otherwise, SLURM will forcefully kill the running container, which can cause the podman state to become corrupted and risks losing the container. If this happens see the instructions HERE The next section is what to do is that happens:

# Debugging 

Sometimes podman stateSuppose the srun session ends while a container is still running, either due to the time limit or typing exit before the container is stopped

Use podman ps to list all containers

podman ps --all

If it lists a container as running, but you are not able to connect to it, then it means the srun session running the container was terminated before the container was stopped

podman exec -it mycontainer /bin/bash
transport endpoint is not connected

At this point the container the podman state is out of sync - run this command to update it:

podman system migrate

By default, both the file system and networking stack of the container are completely isolated from the surrounding environment. This is inconvenient if we wish to access data in the external . In the next section we'll go over how to share information, as well as how to run RStudio or jupyterlab out of a container

# Using containers

You can mount external folders to the container using the -v flag

any changes will be reflected and vice versa

A few more bells and whistles:
- prevents it from terminating
- adds a signal handler

below is the final command:


Helpful links
-------------
https://broad.service-now.com/kb_view.do?sys_kb_id=8923f956479aa91411484438946d4383

https://broad.service-now.com/kb_view.do?sysparm_article=KB0010821

https://backstage.broadinstitute.org/docs/default/component/disco-docs/
