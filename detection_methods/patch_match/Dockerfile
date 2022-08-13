ARG cuda_version=9.0
ARG cudnn_version=7
FROM nvidia/cuda:${cuda_version}-cudnn${cudnn_version}-devel

# Install system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
      bzip2 \
      zip unzip \
      g++ \
      graphviz \
      libgl1-mesa-glx \
      libhdf5-dev \
      openmpi-bin \
      libeigen3-dev \
      libflann-dev \
      libtbb-dev \
      libglu1-mesa-dev freeglut3-dev mesa-common-dev \
      uuid-dev \
      libtbb2 libtbb-dev libjpeg-dev libpng-dev libtiff-dev libjasper-dev libdc1394-22-dev\
      cmake git libgtk2.0-dev pkg-config libavcodec-dev libavformat-dev libswscale-dev \
     libboost-all-dev \
      wget \
      libvtkgdcm2-dev \
      libgdcm-tools\
      libvtkgdcm-tools\
      libvtkgdcm-cil \
      libhighgui-dev\
      libatlas-base-dev \
      libfftw3-dev\
      libblas-dev\
      libomp-dev \
      libinsighttoolkit4-dev &&\
      apt-get install -y --allow-downgrades --allow-change-held-packages &&\
    rm -rf /var/lib/apt/lists/*

# Create folder framework

RUN wget --quiet https://sourceforge.net/projects/opencvlibrary/files/opencv-unix/2.4.13/opencv-2.4.13.zip/download \
          -O opencv-2.4.13.zip && \
          unzip -q opencv-2.4.13.zip && \
          rm  opencv-2.4.13.zip

RUN mkdir -p opencv-2.4.13/release  && \
    cd  opencv-2.4.13/release && \
    cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=/usr/local .. -DWITH_CUDA=OFF && \
    make -j &&  make install && \
    cd ../../ && rm opencv-2.4.13 -fr

COPY forgery-detection-visualization /src/forgery-detection-visualization
WORKDIR /src
RUN cd forgery-detection-visualization && \
    mkdir -p build && cd build && \
    cmake .. &&\
    make -j &&\
    cp bin/forgery_detection /src/forgery_detection -r && \
    cd / && rm forgery-detection -fr 

WORKDIR /src
ENTRYPOINT ["./forgery_detection"]
