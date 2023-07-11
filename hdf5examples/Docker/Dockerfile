FROM fedora:30

ENV LD_LIBRARY_PATH="/usr/lib64/mpich/lib:${LD_LIBRARY_PATH}"
ENV PATH="/usr/lib64/mpich/bin:${PATH}"

RUN dnf install -y hdf5-mpich-devel git cmake rsync git gcc-c++ cmake mpich-devel make gsl-devel eigen3-devel

#RUN dnf install -y make root-mathmore root-mathcore root-minuit2
#RUN dnf install -y root-*
WORKDIR /code
RUN git clone https://github.com/BlueBrain/HighFive.git
RUN git clone https://github.com/diatomic/diy.git


#WORKDIR /data
#COPY sbn_gpvm_example.tar.gz sbn_gpvm_example.tar.gz
#RUN tar xzf sbn_gpvm_example.tar.gz
#COPY small.xml /data/small.xml
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

WORKDIR /
COPY root-6.14.06.tar.gz root-6.14.06.tar.gz
RUN tar xzf root-6.14.06.tar.gz

WORKDIR /code
COPY Prototype.zip Prototype.zip
RUN dnf install -y unzip
RUN dnf install -y freetype
RUN unzip Prototype.zip
#RUN cmake . -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=/usr/local -DDIY_INCLUDE_DIRS=$PWD/diy/include -DHIGHFIVE_DIR=$PWD/HighFive -DHDF5_C_LIBRARIES=/usr/lib64/mpich/lib/libhdf5.so -DHDF5_INCLUDE_DIRS=/usr/include/mpich-x86_64 -DROOTSYS=/root-6.14.06/local -DARCH=knl
RUN cmake . -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=/usr/local -DDIY_INCLUDE_DIRS=$PWD/diy/include -DHIGHFIVE_DIR=$PWD/HighFive -DHDF5_C_LIBRARIES=/usr/lib64/mpich/lib/libhdf5.so -DHDF5_INCLUDE_DIRS=/usr/include/mpich-x86_64 -DROOTSYS=/root-6.14.06/local -DARCH=haswell
RUN make
RUN make install
COPY file /files
ENV LD_LIBRARY_PATH="/root-6.14.06/local/lib:${LD_LIBRARY_PATH}"
