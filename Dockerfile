FROM debian:testing-slim as build

LABEL maintainer="pascal.thibaudeau@cea.fr"
LABEL version="0.2"
LABEL description="Sparselizard on docker. See https://github.com/pthibaud/dockerlizard"

# Disable prompt during packages installation
ARG DEBIAN_FRONTEND=noninteractive

# Update ubuntu software repository
# Install additional repository
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get -y dist-upgrade && \
    apt-get install -y apt-utils build-essential cmake  \
    libopenblas-dev \
    libopenmpi-dev \
    libmumps-dev \
    libmetis-dev \
    petsc-dev \
    slepc-dev \
    libgmsh-dev \
    libomp-dev \
    ninja-build

# Clean the installation
RUN apt-get clean


WORKDIR /sparselizard
COPY CMakeLists.txt .
COPY src src
COPY cmake cmake
COPY simulations simulations
RUN mkdir build
# Prepare compilation environment
RUN cmake -G Ninja -S . -B build
RUN cmake --build build -- -j 8

RUN mkdir -p /usr/local
RUN cmake --install build -v --prefix /usr/local

# Install the library and headers; Put the env variables
# RUN make install
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
