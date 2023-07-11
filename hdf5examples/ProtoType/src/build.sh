cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=/code/src -DEIGEN3_INCLUDE_DIRS=/code/eigen -DDIY_INCLUDE_DIRS=/code/diy/include -DHighFive_INCLUDE_DIR=/code/HighFive/include -DHDF5_C_LIBRARIES=/usr/lib64/mpich/lib/libhdf5.so -DHDF5_INCLUDE_DIRS=/usr/include/mpich-x86_64 -DROOTSYS=/root-6.14.06/local -DMFA_INCLUDE_DIR=/code/mfa/include -DLBFGSB_INCLUDE_DIR=/code/lbfgsb_cpp/include/ -DLBFGSB_LIBRARIES=/code/lbfgsb_cpp/lib/liblbfgsb_cpp.so
make
export LD_LIBRARY_PATH="/root-6.14.06/local/lib:${LD_LIBRARY_PATH}"
make install

