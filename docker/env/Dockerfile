FROM julia:latest

USER root
#add config files
RUN apt-get install unzip
COPY environment.yml Project.toml Manifest.toml docker/env/startup.jl /tmp/
RUN chown $UID:$GID /tmp/environment.yml /tmp/Project.toml /tmp/Manifest.toml /tmp/startup.jl
#COPY postBuild /usr/local/bin/postBuild.sh
#RUN chown $UID:$GID /usr/local/bin/postBuild.sh && \
    #chmod u+x /usr/local/bin/postBuild.sh
COPY docker/entrypoint.sh /usr/local/bin/
RUN chown $UID:$GID /usr/local/bin/entrypoint.sh && \
    chmod u+x /usr/local/bin/entrypoint.sh

USER $USER
#setup startup file
RUN mkdir -p $HOME/.julia/config
RUN cp /tmp/startup.jl $HOME/.julia/config/
# create a project directory inside user home
ENV PROJECT_DIR $HOME/app
RUN mkdir $PROJECT_DIR
RUN cp /tmp/environment.yml /tmp/Manifest.toml /tmp/Project.toml $PROJECT_DIR
WORKDIR $PROJECT_DIR

# build the conda environment
ENV ENV_PREFIX julia-workspace
RUN conda update --name base --channel defaults conda && \
    mamba env create --file $PROJECT_DIR/environment.yml --force && \
    conda clean --all --yes

#build the julia environment TODO find out why RCall does not install properly
ENV R_HOME $CONDA_DIR/envs/$ENV_PREFIX/lib/R
RUN julia --project=$PROJECT_DIR/Project.toml --startup-file=no -e "using Pkg;Pkg.instantiate();Pkg.update()"
RUN julia --project=$PROJECT_DIR/Project.toml --startup-file=no -e "using Pkg;Pkg.precompile()"
# run the postBuild script to install any JupyterLab extensions
#RUN conda activate $ENV_PREFIX && \
 #   /usr/local/bin/postBuild.sh && \
  #  conda deactivate

#setup entrypoint
ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]

# default command will launch JupyterLab server for development
CMD julia --project=$PROJECT_DIR/Project.toml
