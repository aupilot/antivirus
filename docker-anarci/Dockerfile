# syntax=docker/dockerfile:1
FROM ubuntu:20.04

RUN apt-get update && \
    apt-get install -y git python3 python3-pip hmmer && \
    pip install biopython
#    apt install -y curl build-essential git && \

WORKDIR /opt
RUN git clone https://github.com/oxpig/ANARCI.git && \
    cd /opt/ANARCI && \
    python3 setup.py install

#WORKDIR /root
#RUN anarci -h
#RUN /bin/bash


