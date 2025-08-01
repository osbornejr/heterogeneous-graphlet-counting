FROM julia:latest

USER root
#add config files
RUN apt-get update;apt-get install -y unzip ssh binutils locate
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

# build the conda environment (turning off mamba for now)
ENV ENV_PREFIX julia-workspace
#RUN conda update --name base --channel defaults conda=23.3.1
RUN mamba env create --file $PROJECT_DIR/environment.yml && \
    conda clean --all --yes


#build the julia environment TODO find out why RCall does not install properly
ENV R_HOME $CONDA_DIR/envs/$ENV_PREFIX/lib/R
ENV R_LIBS_USER $CONDA_DIR/envs/$ENV_PREFIX/lib/
#ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:$HOME/conda/envs/$ENV_PREFIX/lib
#ENV LD_PRELOAD /usr/lib/x86_64-linux-gnu/libstdc++.so.6
## WORKAROUND AS PER https://discourse.julialang.org/t/glibcxx-3-4-26-not-found-in-rcall/29113/10: need to replace the libstdc++ library from the system (ie what conda uses to determine all other packages including the R dependencies)
RUN cp $HOME/conda/envs/$ENV_PREFIX/lib/libstdc++.so.6 $HOME/julia-${JULIA_VERSION}/lib/julia/
#set up MPI variables
#ENV JULIA_MPI_BINARY system
#ENV JULIA_MPI_PATH $CONDA_DIR/envs/$ENV_PREFIX
#RUN julia --project=$PROJECT_DIR/Project.toml --startup-file=no -e "using Pkg;Pkg.instantiate();Pkg.update()"
#RUN julia --projec#t=$PROJECT_DIR/Project.toml --startup-file=no -e "using Pkg; for pkg in collect(keys(Pkg.installed()));Base.compilecache(Base.identify_package(pkg));end"
#set up MPI.jl to use openMPI 
#RUN julia --project=$PROJECT_DIR/Project.toml --startup-file=no -e "using Pkg;Pkg.build(\"MPI\");using MPI;MPI.install_mpiexecjl()"
#Or without MPI-- TODO fix it so that can instantiate during build... currently we must run instantiate everytime we install on a fresh machine- it does not work here currently because the dev packages will not install before they are linked to in a container. Could solve maybe by registering them?
COPY packages/ $PROJECT_DIR/packages/
RUN julia --project=$PROJECT_DIR/Project.toml --startup-file=no -e "using Pkg;Pkg.instantiate();Pkg.update()"
#add julia bin to path
ENV  PATH $PATH:$HOME/.julia/bin

# run the postBuild script to install any JupyterLab extensions
#RUN conda activate $ENV_PREFIX && \
 #   /usr/local/bin/postBuild.sh && \
  #  conda deactivate

##set up env variable to dynamically change startupfile to be loaded by julia with -i flag
ENV startup docker/env/startup.jl
#setup entrypoint
ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]

# default command will launch julia
CMD ["sh", "-c", "julia --project=$PROJECT_DIR/Project.toml -i ${startup}"]
