# Docker for CONCOCT (http://github.com/BinPro/CONCOCT) v1.1.0
# VERSION 1.1.0
#
# This docker creates and sets up a Debian environment with all
# dependencies for CONCOCT v1.1.0 installed.
#
# To login to the docker with a shared directory from the host do:
#
# docker run -v /my/host/shared/directory:/my/docker/location -i -t alneberg/concoct_1.1.0 /bin/bash
#

FROM python:3.12-slim-bookworm
COPY . /opt/CONCOCT

# Get basic packages needed
RUN apt-get update -qq && \
    apt-get install -qq --no-install-recommends wget ca-certificates build-essential libgsl-dev git zip unzip bedtools samtools && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/BinPro/integration_test_data/archive/v1.0.tar.gz && \
    mkdir /opt/CONCOCT/tests/test_data/integration_test_data && \
    tar -xzf v1.0.tar.gz -C /opt/CONCOCT/tests/test_data/integration_test_data --strip-components=1 && \
    rm v1.0.tar.gz

# Install python dependencies and install CONCOCT 1.1.0
RUN cd /opt/CONCOCT && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir .

RUN cd /opt/CONCOCT && \
    pytest
