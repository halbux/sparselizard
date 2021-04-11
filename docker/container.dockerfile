# syntax=docker/dockerfile:1
FROM ubuntu:20.04
RUN apt-get update
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
RUN apt-get update && apt-get install -y \
    build-essential \
    gpg \
    software-properties-common \
    ca-certificates \
    wget \
    curl \
    python3-distutils \
    python3 \
    git \
    gfortran


RUN apt purge --auto-remove cmake -y
RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null
RUN apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
RUN apt update
RUN apt install kitware-archive-keyring
RUN rm /etc/apt/trusted.gpg.d/kitware.gpg
# GCC/G++ 11
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys DE19EB17684BA42D
RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y

RUN apt update

RUN apt install cmake g++-11 gcc-11 gfortran-11 -y

RUN rm -f /usr/bin/gcc /usr/bin/g++ /usr/bin/gfortran
RUN ln -s /usr/bin/gcc-11 /usr/bin/gcc
RUN ln -s /usr/bin/g++-11 /usr/bin/g++
RUN ln -s /usr/bin/gfortran-11 /usr/bin/gfortran

# pass only the dependency script to cache the docker image layer across multiple builds
COPY install_external_libs/install_petsc_auto.sh /root/
RUN chmod 777 ~/install_petsc_auto.sh
RUN bash ~/install_petsc_auto.sh

# pass only the dependency script to cache the docker image layer across multiple builds
COPY install_external_libs/install_gmsh.sh /root/
RUN apt install gmsh -y
RUN bash ~/install_gmsh.sh


COPY . /app
WORKDIR /app

CMD bash mnrun.sh 5 55