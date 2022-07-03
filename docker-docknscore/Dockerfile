FROM nvidia/cuda:11.4.1-cudnn8-devel-ubuntu20.04

RUN apt-key del 7fa2af80 && \
    apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/3bf863cc.pub && \
    apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/machine-learning/repos/ubuntu2004/x86_64/7fa2af80.pub

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Australia/Brisbane
RUN apt-get update && apt-get install -y tzdata
RUN apt-get install -y --no-install-recommends git python3-dev python3-pip cmake libssl-dev libncurses5-dev libncursesw5-dev doxygen swig software-properties-common wget mc

RUN python3 -m pip install --upgrade pip && \
    pip install numpy

RUN apt-get install -y cuda-samples-11-4 lz4

# ############# FFTW ################
WORKDIR /tmp
ENV FFTW_TARGET fftw-3.3.10
RUN wget -P /tmp http://www.fftw.org/${FFTW_TARGET}.tar.gz && \
    tar xzf /tmp/${FFTW_TARGET}.tar.gz -C /tmp && \
    rm -rf /tmp/${FFTW_TARGET}.tar.gz && \
    cd /tmp/${FFTW_TARGET} && \
    ./configure --enable-float --enable-avx2 && \
    make -j$(nproc) && \
    make install

################ DLA-Ranker/kir ########################
#RUN pip install torch torchvision --extra-index-url https://download.pytorch.org/whl/cu113
RUN pip install matplotlib lz4 tensorflow-gpu==2.9.0
RUN pip install prody
RUN pip install scikit-learn==0.22.1 pandas

WORKDIR /opt
RUN git clone https://github.com/aupilot/DLA-Ranker.git && \
    rm -Rf /opt/DLA-Ranker/Models/CCD4PPI && \
    rm -Rf /opt/DLA-Ranker/Models/BM5

# we need these to run naccess crap
RUN apt-get install -y csh libgfortran5

############# MEGADOCK KIR #############
# important: use branch "kir"!
# mega-post.sh file in repo is invalid, we copy it over later
ENV MEGADOCK_WORK_DIR /opt/megadock-kir
RUN git clone -b kir https://github.com/aupilot/MEGADOCK.git ${MEGADOCK_WORK_DIR}
RUN cp -f ${MEGADOCK_WORK_DIR}/Dockerfiles/gpu/Makefile ${MEGADOCK_WORK_DIR}/Makefile
RUN cd ${MEGADOCK_WORK_DIR} && \
    make -j$(nproc)

WORKDIR ${MEGADOCK_WORK_DIR}
RUN ln -s ${MEGADOCK_WORK_DIR}/megadock-gpu /usr/local/bin/megadock-gpu && \
    ln -s ${MEGADOCK_WORK_DIR}/decoygen /usr/local/bin/decoygen && \
    ln -s ${MEGADOCK_WORK_DIR}/block /usr/local/bin/block && \
    ln -s ${MEGADOCK_WORK_DIR}/ppiscore /usr/local/bin/ppiscore && \
    ln -s ${MEGADOCK_WORK_DIR}/mega-post.sh /usr/local/bin/mega-post.sh

COPY mega-post.sh /opt/megadock-kir/mega-post.sh

########### Open Babel ###########
# open babael official repo has broken
# we have a fixed copy, so go the long way
RUN apt-get install -y libeigen3-dev libxml2
WORKDIR /tmp
RUN git clone https://github.com/aupilot/openbabel.git
RUN mkdir /tmp/openbabel/build
WORKDIR /tmp/openbabel/build
RUN cmake -DBUILD_GUI=NO .. && make -j8 && make install

############## Vina ##############
#ADD https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.3/vina_1.2.3_linux_x86_64 /usr/local/bin/vina
ADD https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.1.2-boost-new/vina_1.1.2-boost-new_linux_x86_64 /usr/local/bin/vina
RUN chmod +x /usr/local/bin/vina
# we don't use python integration
#RUN pip install vina

########### OnionNet #############
# Do we need that specific version??   scikit-learn==0.22.1
RUN pip install rdkit biopandas joblib mdtraj
RUN git clone -b master https://github.com/aupilot/onionnet /opt/onionnet
#ADD https://drive.google.com/uc?export=download&id=14QmxihpshNpxwev2Cvnw2Qadspr1WLbF&confirm=Y /opt/onionnet/models_new/CNN_final_model_weights.h5
ADD http://dl.dropboxusercontent.com/s/w6kl81cr84ayniz/CNN_final_model_weights.h5?dl=1 /opt/onionnet/models_new/CNN_final_model_weights.h5
#RUN mkdir -p /opt/onionnet/models_new
#COPY CNN_final_model_weights.h5v /opt/onionnet/models_new/CNN_final_model_weights.h5

RUN mkdir -p /opt/var
COPY dock-n-score.py /opt/dock-n-score.py
COPY block3.py /opt/block3.py
COPY run_obabel.sh /opt/run_obabel.sh
#COPY onion_score.sh /opt/onionnet/kir01/

RUN pip cache purge && \
    apt-get clean

# TODO: we might need to install the original megadock too?