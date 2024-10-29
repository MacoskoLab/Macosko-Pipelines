Intro to containers
-------------------

Let's start a new srun session, this time giving us more time and memory:

```
srun --partition=hpcx_macosko --time 1-00:00:00 --mem 32G --pty bash
```

Finally, it's time to make a container. This command creates a new interactive container:

```
podman create -it --name mycontainer --init myimage
```

You can list all podman containers with this command:

```
podman ps -a
```

You should see a new container with the status "Created". Now let's start and attach to it:

```
podman start -a mycontainer
```

You are now running a bash terminal inside the container! By default, you will be root with home directory /root. If you take a look around with `ls` you'll see that there is a whole new filesystem in here. Exit anytime:

```
exit
```

The podman container automatically stops when the main process exits, regardless if there are background jobs running. If you started the container without `-a` you can manually stop it using `podman stop -t 0 mycontainer`.

In addition to `/local/podman/containers-user-$UID`, container data is stored in `/tmp/containers-user-$UID` and running containers are stored in `/tmp/podman-run-$UID`. Therefore, like images, containers are stored in node-specific directories and will not be available on other nodes.


We are done using this container as an example - remove it like so

```
podman rm mycontainer
```

While it's possible to treat a container as a compute instance - where you download your data and install packages into its filesystem, starting/stopping it each time you need it - this kind of stateful "pet" container is generally discouraged. Containers should be treated as immutable for the most part, so it's better to store your files and environment either as part of the image or on the host filesystem. Then, each time you want to start a session, you creating a new container from the base image.

This commands creates a new container, allows you to run code inside of it, then removes itself when completed:

```
podman run -it --rm --name mycontainer --init myimage
```

By default, both the file system and networking stack of the container are completely isolated from the surrounding host environment.

You can bind mount external folders to the container using the -v flag:

```
podman run -it --rm --name mycontainer -v /broad/macosko/$USER:/broad/macosko/$USER --init myimage
```

You can also map ports from the container to the host:

```
podman run -it --rm --name mycontainer -p 8787:8787 --init myimage
```

See the Quick Start guide for a set of bash commands that automates the process of finding ports and starting containers

If you've made a lot of changes to the container and want to use it as the base for next time, you can commit it to an image:

```
podman commit myimage newimage
```

Helpful links
-------------
https://broad.service-now.com/kb_view.do?sys_kb_id=8923f956479aa91411484438946d4383

https://broad.service-now.com/kb_view.do?sysparm_article=KB0010821

https://backstage.broadinstitute.org/docs/default/component/disco-docs/using-containers/
