#!/bin/bash --login
set -e
conda activate $ENV_PREFIX
#setup workspace
#jupyter lab workspaces import workspace.json
#set up password for notebook
#jupyter notebook --generate-config
#echo "c.NotebookApp.token='jupyter-julia'">> ~/.jupyter/jupyter_notebook_config.py
exec "$@"
