FROM ubuntu:jammy

RUN apt-get update && \
    apt-get install -y \
        cmake \
        gcc \
        g++ \
        libx11-dev \
        libglu1-mesa-dev \
        xorg-dev \
        libglfw3-dev \
    && rm -rf /var/lib/apt/lists/*


#--- Install MPICH

ARG mpich_version
RUN apt-get update && \
    apt-get install -y \
        mpich${mpich_version+=${mpich_version:?}} \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    apt-get install -y \
        git \
    && rm -rf /var/lib/apt/lists/*
