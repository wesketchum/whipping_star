FROM iamholger/hpcbase:0.1

WORKDIR /build
RUN git clone https://github.com/steven-varga/h5cpp.git
RUN git clone https://github.com/diatomic/diy.git
RUN git clone https://github.com/BlueBrain/HighFive.git

WORKDIR /build/h5cpp
COPY sbn-diyh5cpp-0.0.1.tar.gz sbn-diyh5cpp-0.0.1.tar.gz
RUN tar xzf sbn-diyh5cpp-0.0.1.tar.gz
RUN cmake sbn-diyh5cpp-0.0.1 -DCMAKE_INSTALL_PREFIX=/usr/local -DDIY_INCLUDE_DIRS=/build/diy/include -DH5CPP_DIR=/build/h5cpp -DHDF5_C_LIBRARIES=/usr/local/lib/libhdf5.so -DHDF5_INCLUDE_DIRS=/usr/local/include
RUN make
RUN make install

WORKDIR /build/HighFive
COPY sbn-diyhighfive-0.0.1.tar.gz sbn-diyhighfive-0.0.1.tar.gz
RUN tar xzf sbn-diyhighfive-0.0.1.tar.gz
RUN cmake sbn-diyhighfive-0.0.1 -DCMAKE_INSTALL_PREFIX=/usr/local -DDIY_INCLUDE_DIRS=/build/diy/include -DHIGHFIVE_DIR=/build/HighFive -DHDF5_C_LIBRARIES=/usr/local/lib/libhdf5.so -DHDF5_INCLUDE_DIRS=/usr/local/include
RUN make
RUN make install
