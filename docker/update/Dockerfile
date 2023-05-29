FROM osbornejr/julia-env:latest
USER root
#add config files
COPY environment.yml Project.toml Manifest.toml docker/env/startup.jl /tmp/
RUN chown $UID:$GID /tmp/environment.yml /tmp/Project.toml /tmp/Manifest.toml /tmp/startup.jl
#COPY postBuild /usr/local/bin/postBuild.sh
#RUN chown $UID:$GID /usr/local/bin/postBuild.sh && \
    #chmod u+x /usr/local/bin/postBuild.sh
COPY docker/entrypoint.sh /usr/local/bin/
RUN chown $UID:$GID /usr/local/bin/entrypoint.sh && \
    chmod u+x /usr/local/bin/entrypoint.sh

#update environments
USER $USER
#setup startup file
RUN mkdir -p $HOME/.julia/config
RUN cp /tmp/startup.jl $HOME/.julia/config/
RUN cp /tmp/environment.yml /tmp/Manifest.toml /tmp/Project.toml $PROJECT_DIR
WORKDIR $PROJECT_DIR

# build the conda environment
ENV ENV_PREFIX julia-workspace
#RUN conda update --name base --channel defaults conda
#RUN conda update -n base -c defaults conda
RUN mamba env update --file $PROJECT_DIR/environment.yml && \
    conda clean --all --yes

#build the julia environment TODO find out why RCall does not install properly
ENV R_HOME $CONDA_DIR/envs/$ENV_PREFIX/lib/R
ENV R_LIBS_USER $CONDA_DIR/envs/$ENV_PREFIX/lib/


#RUN julia --project=$PROJECT_DIR/Project.toml --startup-file=no -e "using Pkg;Pkg.instantiate();Pkg.update()"
##TODO FIX ISSUE WHERE DEV PACKAGES (IN Project.toml) CANT BE PRECOMPILED HERE
#RUN julia --project=$PROJECT_DIR/Project.toml -e "using Pkg; for pkg in collect(keys(Pkg.installed()));Base.compilecache(Base.identify_package(pkg));end" 
#RUN julia --project=$PROJECT_DIR/Project.toml -e "using Pkg;using CSV"
# run the postBuild script to install any JupyterLab extensions
#RUN conda activate $ENV_PREFIX && \
 #   /usr/local/bin/postBuild.sh && \
  #  conda deactivate

ENV startup docker/env/startup.jl
#setup entrypoint
ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]

# default command will launch julia
CMD ["sh", "-c", "julia --project=$PROJECT_DIR/Project.toml -i ${startup}"]

