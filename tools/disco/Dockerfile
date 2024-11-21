FROM pipeline-image:latest

ARG USER=REQUIRED
RUN [ "$USER" != "REQUIRED" ]

### Configure RStudio ###
# Set password
RUN echo "root:$USER" | chpasswd
# Enable root to run RStudio
RUN echo "auth-minimum-user-id=0" >> /etc/rstudio/rserver.conf
# Set the working directory and package directory to /broad/macosko/$USER
RUN echo "wd <- '/broad/macosko/$USER'\nif (dir.exists(wd)) {\n  setwd(wd)\n  pd <- file.path(wd, 'R-packages')\n  if (!dir.exists(pd)) {dir.create(pd)}\n  .libPaths(c(pd, .libPaths()))\n  if (.libPaths()[1] != pd) {print('WARNING: Could not set user package directory')}\n} else {\n  print('WARNING: User directory not found')\n}\nsuppressWarnings(rm(wd,pd))" > /root/.Rprofile
# (miscellaneous config)
RUN echo '{\n  "sync_files_pane_working_dir": true\n}' > /etc/rstudio/rstudio-prefs.json

### Configure JupyterLab ###
# Set password
RUN expect -c 'spawn bash -lc "micromamba run jupyter lab password"; expect "Enter password:"; send "$env(USER)\r"; expect "Verify password:"; send "$env(USER)\r"; expect eof'
# Set the mamba root to /broad/macosko/$USER/micromamba if available
RUN sed -i "4i\  [ ! -d '/broad/macosko/$USER' ] && echo 'WARNING: User directory not found' || bash -ic 'micromamba shell init -s bash -r /broad/macosko/$USER/micromamba' &> /dev/null" /root/.profile
# (miscellaneous config)
RUN /bin/bash -lc "micromamba run jupyter lab --generate-config && echo -e 'c.ServerApp.terminals_enabled = True\nc.FileContentsManager.delete_to_trash = False' >> /root/.jupyter/jupyter_lab_config.py"

### Create bash functions ###
RUN printf "\nrstudio() { \n\
    pgrep -x rserver > /dev/null && echo 'rserver is already running' && return 1 \n\
    rstudio-server start \n\
    exec bash \n}\n" >> /root/.bashrc
RUN printf "\njupyterlab() { \n\
    [ ! -d \$MAMBA_ROOT_PREFIX ] && echo 'ERROR: mamba root could not be found' && return 1 \n\ 
    pgrep -x jupyter-lab > /dev/null && echo 'jupyter-lab is already running' && return 1 \n\
    micromamba list | grep -q jupyterlab || micromamba install -y jupyterlab python=3.11.2 \n\
    tmux new -s jupyterlab -d bash -ic \"micromamba run -n ${1:-base} jupyter lab --allow-root --ip='*' --port='8787' --NotebookApp.token='' --no-browser\" \n\
    exec bash \n}\n" >> /root/.bashrc

ENTRYPOINT ["/bin/bash", "-lc"]
CMD ["/bin/bash", "-i"]
