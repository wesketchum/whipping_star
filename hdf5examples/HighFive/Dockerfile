FROM iamholger/hpcbase:0.1

WORKDIR /build
COPY sbn-diyh5cpp-0.0.1.tar.gz sbn-diyh5cpp-0.0.1.tar.gz
RUN tar xzf sbn-diyh5cpp-0.0.1.tar.gz
RUN git clone https://github.com/steven-varga/h5cpp.git
RUN git clone https://github.com/diatomic/diy.git
RUN cmake sbn-diyh5cpp-0.0.1 -DCMAKE_INSTALL_PREFIX=/usr/local -DDIY_INCLUDE_DIRS=$PWD/diy/include -DH5CPP_DIR=$PWD/h5cpp -DHDF5_C_LIBRARIES=/usr/local/lib/libhdf5.so -DHDF5_INCLUDE_DIRS=/usr/local/include
RUN make
RUN make install
