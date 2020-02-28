FROM ubuntu:latest

RUN apt-get update && apt-get -y upgrade
RUN apt-get -y install wget curl unzip build-essential zlib1g-dev git python3 python3-pip openjdk-8-jre libbz2-dev libssl-dev cython

RUN pip3 install cerebra