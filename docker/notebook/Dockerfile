#Dockerfile created following this guide: https://towardsdatascience.com/conda-pip-and-docker-ftw-d64fe638dc45
FROM ubuntu:18.04

SHELL [ "/bin/bash", "--login", "-c" ]
# Create a non-root user
ARG username=osbornejr
ARG uid=1000
ARG gid=100
ENV USER $username
ENV UID $uid
ENV GID $gid
ENV HOME /home/$USER
RUN adduser --disabled-password \
    --gecos "Non-root user" \
    --uid $UID \
    --gid $GID \
    --home $HOME \
    $USER

#add config files
COPY environment.yml Project.toml Manifest.toml /tmp/
RUN chown $UID:$GID /tmp/environment.yml /tmp/Project.toml /tmp/Manifest.toml
#COPY postBuild /usr/local/bin/postBuild.sh
#RUN chown $UID:$GID /usr/local/bin/postBuild.sh && \
    #chmod u+x /usr/local/bin/postBuild.sh
COPY docker/entrypoint.sh /usr/local/bin/
RUN chown $UID:$GID /usr/local/bin/entrypoint.sh && \
    chmod u+x /usr/local/bin/entrypoint.sh

#install dependencies
RUN apt-get -qq update && apt-get -qq -y install curl wget bzip2 

#switch to non-root user
USER $USER
ENV CONDA_DIR $HOME/conda
##install miniconda (adapted from official miniconda container)
RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o $HOME/miniconda.sh \
    && bash $HOME/miniconda.sh -bfp $CONDA_DIR \
    && rm -rf $HOME/miniconda.sh
# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile
# make conda activate command available from /bin/bash --interative shells
RUN conda init bash \
    && conda update conda


USER root
# Julia dependencies
# install Julia packages in /opt/julia instead of $HOME
#ENV JULIA_DEPOT_PATH=/opt/julia
#ENV JULIA_PKGDIR=/opt/julia
ENV JULIA_VERSION=1.4.1

WORKDIR /tmp

# hadolint ignore=SC2046
RUN mkdir "$HOME/julia-${JULIA_VERSION}" && \
    wget -q https://julialang-s3.julialang.org/bin/linux/x64/$(echo "${JULIA_VERSION}" | cut -d. -f 1,2)"/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" && \
    echo "fd6d8cadaed678174c3caefb92207a3b0e8da9f926af6703fb4d1e4e4f50610a *julia-${JULIA_VERSION}-linux-x86_64.tar.gz" | sha256sum -c - && \
    tar xzf "julia-${JULIA_VERSION}-linux-x86_64.tar.gz" -C "$HOME/julia-${JULIA_VERSION}" --strip-components=1 && \
    rm "/tmp/julia-${JULIA_VERSION}-linux-x86_64.tar.gz"
RUN ln -fs $HOME/julia-*/bin/julia /usr/local/bin/julia

## Show Julia where conda libraries are \
#RUN mkdir /etc/julia && \
#    echo "push!(Libdl.DL_LOAD_PATH, \"$CONDA_DIR/lib\")" >> /etc/julia/juliarc.jl && \
#    # Create JULIA_PKGDIR \
#    mkdir "${JULIA_PKGDIR}" && \
#    chown "${NB_USER}" "${JULIA_PKGDIR}"
#   # fix-permissions "${JULIA_PKGDIR}"

#remove dependencies
RUN apt-get -qq -y remove curl bzip2\
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

USER $USER
##Setting up conda
# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH
# make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile
# make conda activate command available from /bin/bash --interative shells
RUN conda init bash

# create a project directory inside user home
ENV PROJECT_DIR $HOME/app
RUN mkdir $PROJECT_DIR
RUN cp /tmp/environment.yml /tmp/Manifest.toml /tmp/Project.toml $PROJECT_DIR
WORKDIR $PROJECT_DIR

# build the conda environment
ENV ENV_PREFIX julia-workspace
RUN conda update --name base --channel defaults conda && \
    conda env create --file $PROJECT_DIR/environment.yml --force && \
    conda clean --all --yes

#build the julia environment TODO find out why RCall does not install properly
ENV R_HOME $CONDA_DIR/envs/$ENV_PREFIX/lib/R
RUN julia --project=$PROJECT_DIR/Project.toml -e "using Pkg;Pkg.instantiate();Pkg.precompile()"
# run the postBuild script to install any JupyterLab extensions
#RUN conda activate $ENV_PREFIX && \
 #   /usr/local/bin/postBuild.sh && \
  #  conda deactivate

#setup entrypoint
ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]

# default command will launch JupyterLab server for development
CMD [ "jupyter", "lab", "--no-browser", "--ip", "0.0.0.0" ]
