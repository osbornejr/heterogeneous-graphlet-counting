#!/bin/bash --login
set -e
conda activate $ENV_PREFIX
#set number of threads for julia to use.
export JULIA_NUM_THREADS=$(nproc)
#setup workspace
#jupyter lab workspaces import workspace.json
#set up password for notebook
#jupyter notebook --generate-config
#echo "c.NotebookApp.token='jupyter-julia'">> ~/.jupyter/jupyter_notebook_config.py
exec "$@"
