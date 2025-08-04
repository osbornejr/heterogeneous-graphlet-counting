FROM conda-base:latest

USER root
RUN apt-get update && apt-get install -y \
	wget 

# Julia dependencies
# install Julia packages in /opt/julia instead of $HOME
#ENV JULIA_DEPOT_PATH=/opt/julia
#ENV JULIA_PKGDIR=/opt/julia
ENV JULIA_VERSION=1.11.6
WORKDIR /tmp
ARG TARGETARCH
RUN case "$TARGETARCH" in \
    "amd64")  export JULIA_ARCH="x86_64";; \
    "arm64")  export JULIA_ARCH="aarch64";; \
    *) echo "Unsupported arch"; exit 1;; \
    esac && \
    mkdir "$HOME/julia-${JULIA_VERSION}" && \
    wget -q https://julialang-s3.julialang.org/bin/linux/${JULIA_ARCH}/$(echo "${JULIA_VERSION}" | cut -d. -f 1,2)"/julia-${JULIA_VERSION}-linux-${JULIA_ARCH}.tar.gz" && \
   # echo "f190c938dd6fed97021953240523c9db448ec0a6760b574afd4e9924ab5615f1 *julia-${JULIA_VERSION}-linux-x86_64.tar.gz" | sha256sum -c - && \
    tar xzf "julia-${JULIA_VERSION}-linux-${JULIA_ARCH}.tar.gz" -C "$HOME/julia-${JULIA_VERSION}" --strip-components=1 && \
    rm "/tmp/julia-${JULIA_VERSION}-linux-${JULIA_ARCH}.tar.gz"
RUN ln -fs $HOME/julia-*/bin/julia /usr/local/bin/julia

USER $USER
CMD ["julia"]



