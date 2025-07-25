FROM ubuntu:20.04

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

#detect target architecture
ARG TARGETARCH
RUN echo "TARGETARCH is $TARGETARCH"

# Set  variable conditionally during build
RUN if [ "$TARGETARCH" = "amd64" ]; then \
        echo "CONDAARCH=x86_64" >> /tmp/env.txt; \
    elif [ "$TARGETARCH" = "arm64" ]; then \
        echo "CONDAARCH=aarch64" >> /tmp/env.txt; \
    else \
        echo "CONDAARCH=x86_64" >> /tmp/env.txt; \
    fi 

#install dependencies
RUN apt-get update && apt-get install -y \
	curl 

#add config files
#COPY environment.yml Project.toml Manifest.toml /tmp/
#RUN chown $UID:$GID /tmp/environment.yml /tmp/Project.toml /tmp/Manifest.toml
#COPY postBuild /usr/local/bin/postBuild.sh
#RUN chown $UID:$GID /usr/local/bin/postBuild.sh && \
    #chmod u+x /usr/local/bin/postBuild.sh
#COPY docker/entrypoint.sh /usr/local/bin/
#RUN chown $UID:$GID /usr/local/bin/entrypoint.sh && \
#    chmod u+x /usr/local/bin/entrypoint.sh

#switch to non-root user
USER $USER
##install miniconda (adapted from official miniconda container)
ENV CONDA_DIR $HOME/conda
#RUN . /tmp/env.txt && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-$CONDAARCH.sh -o $HOME/miniconda.sh \
RUN . /tmp/env.txt && curl -sSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -o $HOME/miniconda.sh \
    && bash $HOME/miniconda.sh -bfp $CONDA_DIR \
    && rm -rf $HOME/miniconda.sh
# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile
# make conda activate command available from /bin/bash --interative shells
RUN conda init bash
RUN conda --version
RUN conda update --name base --channel conda-forge conda
RUN conda --version
RUN conda install -c conda-forge mamba

USER root
RUN apt-get clean

USER $USER
	
CMD [ "/bin/bash" ]
