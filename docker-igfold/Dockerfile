FROM nvidia/cuda:11.4.1-cudnn8-devel-ubuntu20.04

RUN apt-key del 7fa2af80 && \
    apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/3bf863cc.pub && \
    apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/machine-learning/repos/ubuntu2004/x86_64/7fa2af80.pub

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Australia/Brisbane
RUN apt-get update && apt-get install -y tzdata
RUN apt-get install -y git python3 python3-pip cmake libssl-dev libncurses5-dev libncursesw5-dev doxygen swig software-properties-common wget mc

WORKDIR /opt
RUN wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1.tar.gz && \
    tar -xf cmake-3.23.1.tar.gz && \
    cd cmake-3.23.1 && \
    ./bootstrap && \
    make -j8 && \
    make install

RUN python3 -m pip install --upgrade pip && \
    pip install  numpy Cython
#protobuf==3.9.2

RUN pip3 install torch torchvision --extra-index-url https://download.pytorch.org/whl/cu113

WORKDIR /opt
RUN git clone https://github.com/openmm/openmm.git && \
    cd openmm && \
    git checkout 7.6.0 && \
    mkdir build && \
    cd build  && \
    cmake .. -DOPENMM_BUILD_CUDA_DOUBLE_PRECISION_TESTS=OFF && \
    make -j8 && \
    make install && \
    make PythonInstall

WORKDIR /opt
RUN git clone https://github.com/openmm/pdbfixer.git  && \
    cd pdbfixer && \
    python3 setup.py install

RUN pip install igfold

ADD https://data.graylab.jhu.edu/IgFold.tar.gz /usr/local/lib/python3.8/dist-packages/igfold/trained_models/

WORKDIR /usr/local/lib/python3.8/dist-packages/igfold/trained_models/
RUN tar -xzf IgFold.tar.gz && \
    rm IgFold.tar.gz

RUN pip cache purge && \
    apt-get clean && \
    rm -Rf /opt/cmake-3.23.1 && \
    rm -Rf /opt/openmm

COPY ig_fold.py /opt/IgFold/ig_fold.py
