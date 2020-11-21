FROM debian:testing

# install fonts and basic programs
RUN apt-get update -q \
        && DEBIAN_FRONTEND=noninteractive \
        apt-get install -qy --no-install-recommends --no-install-suggests \
        g++ \
        make \
        cmake \
        libopenblas-dev \
        libgmsh-dev \
        libmetis-dev \
        libopenmpi-dev \
        libpetsc-real-dev \
        libslepc-real-dev \
        libmumps-dev \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

# source directory shell be mounted here
WORKDIR /workdir
VOLUME /workdir

RUN useradd --uid 1000 --no-create-home --shell /bin/bash builduser
USER builduser

ENTRYPOINT docker/buildscript.sh
