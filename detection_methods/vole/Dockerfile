ARG cuda_version=9.0
ARG cudnn_version=7
FROM nvidia/cuda:${cuda_version}-cudnn${cudnn_version}-devel

# Install system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
      bzip2 \
      zip unzip \
      g++ \
      graphviz \
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


# Download and install opencv-2.4
RUN wget --quiet https://sourceforge.net/projects/opencvlibrary/files/opencv-unix/2.4.13/opencv-2.4.13.zip/download \
          -O opencv-2.4.13.zip && \
          unzip -q opencv-2.4.13.zip && \
          rm  opencv-2.4.13.zip

RUN mkdir -p opencv-2.4.13/release  && \
    cd  opencv-2.4.13/release && \
    cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=/usr/local .. -DWITH_CUDA=OFF && \
    make -j &&  make install && \
    cd ../../ && rm opencv-2.4.13 -fr

# Install vole
ADD ./cmfd_framework/ /src/cmfd_framework
RUN cd /src/cmfd_framework && mkdir -p build && cd build && \
    cmake -DVole_CMFD=ON -DVole_CMFD_Ground_Truth=ON -DVole_Shell=ON .. &&\
    make -j &&\
    cp bin/vole /src/vole &&\
    cd /src && rm -fr cmfd_framework 


WORKDIR /src
ENTRYPOINT ["./vole"]
