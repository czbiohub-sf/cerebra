# run as:
#   docker build -t cerebra .

#start off with a plain Debian
FROM ubuntu:latest

# basic setup stuff
RUN apt-get update
RUN apt-get -y upgrade --fix-missing
RUN apt-get -y install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev git python3 python3-pip

#cerebra proper
COPY . /cerebra
RUN cd /cerebra && pip3 install -r requirements.txt && python3 setup.py install
