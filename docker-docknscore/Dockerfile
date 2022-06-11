FROM nvidia/cuda:11.4.1-cudnn8-devel-ubuntu20.04

RUN apt-key del 7fa2af80 && \
    apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/3bf863cc.pub && \
    apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/machine-learning/repos/ubuntu2004/x86_64/7fa2af80.pub

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Australia/Brisbane
RUN apt-get update && apt-get install -y tzdata
RUN apt-get install -y --no-install-recommends git python3 python3-pip cmake libssl-dev libncurses5-dev libncursesw5-dev doxygen swig software-properties-common wget mc

RUN python3 -m pip install --upgrade pip && \
    pip install numpy


RUN apt-get update -y
RUN apt-get install -y cuda-samples-11-4

# FFTW
ENV FFTW_TARGET fftw-3.3.10

RUN wget -P /tmp http://www.fftw.org/${FFTW_TARGET}.tar.gz && \
    tar xzf /tmp/${FFTW_TARGET}.tar.gz -C /tmp && \
    rm -rf /tmp/${FFTW_TARGET}.tar.gz && \
    cd /tmp/${FFTW_TARGET} && \
    ./configure --enable-float --enable-avx2 && \
    make -j$(nproc) && \
    make install

# MEGADOCK
WORKDIR /opt
RUN git clone https://github.com/aupilot/MEGADOCK.git

ENV MEGADOCK_WORK_DIR /opt/MEGADOCK

RUN cp -f ${MEGADOCK_WORK_DIR}/Dockerfiles/gpu/Makefile ${MEGADOCK_WORK_DIR}/Makefile

RUN cd ${MEGADOCK_WORK_DIR} && \
    make -j$(nproc)


#RUN pip3 install torch torchvision --extra-index-url https://download.pytorch.org/whl/cu113

#RUN pip cache purge && \
#    apt-get clean && \
#    rm -Rf /opt/cmake-3.23.1 && \
#    rm -Rf /opt/openmm

COPY dock-n-score.py /opt/dock-n-score.py
COPY block3.py /opt/block3.py
COPY run_obabel.sh /opt/run_obabel.sh



