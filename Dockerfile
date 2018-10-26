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
RUN apt-get install -qq wget build-essential libgsl0-dev git zip unzip bedtools python-pip

RUN pip install --upgrade pip

COPY . /opt/CONCOCT

# Install python dependencies and fetch and install CONCOCT 0.5.0
RUN cd /opt/CONCOCT;\
    pip install -r requirements.txt;

RUN cd /opt/CONCOCT/;\
    python setup.py install

RUN cd /opt/CONCOCT/;\
    nosetests
