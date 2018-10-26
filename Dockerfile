# Docker for CONCOCT (http://github.com/BinPro/CONCOCT) v0.5.0
# VERSION 0.5.0
#
# This docker creates and sets up an Ubuntu environment with all
# dependencies for CONCOCT v0.5.0 installed.
#
# To login to the docker with a shared directory from the host do:
#
# docker run -v /my/host/shared/directory:/my/docker/location -i -t alneberg/concoct_0.5.0 /bin/bash
#

FROM ubuntu:18.04

# Get basic ubuntu packages needed
RUN apt-get update -qq
RUN apt-get install -qq wget build-essential libgsl0-dev git zip unzip bedtools

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

RUN cd /opt;\
    wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;\
    chmod +x miniconda.sh;\
    ./miniconda.sh -p /opt/miniconda -b

# Use python 3 from conda
ENV PATH="/opt/miniconda/bin:$PATH"

COPY . /opt/CONCOCT

# Install python dependencies and fetch and install CONCOCT 0.5.0
RUN cd /opt/CONCOCT;\
    pip install -r requirements.txt;

RUN cd /opt/CONCOCT/;\
    python setup.py install

RUN cd /opt/CONCOCT/;\
    nosetests
