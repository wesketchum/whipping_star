cmake . -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=/usr/local -DMFA_INCLUDE_DIR=$PWD/mfa/include -DDIY_INCLUDE_DIRS=$PWD/diy/include -DEIGEN3_INCLUDE_DIRS=$PWD/eigen -DLBFGS_INCLUDE_DIRS=$PWD/LBFGSpp -DHighFive_INCLUDE_DIR=$PWD/HighFive/include -DHDF5_C_LIBRARIES=/usr/lib64/mpich/lib/libhdf5.so -DHDF5_INCLUDE_DIRS=/usr/include/mpich-x86_64 -DROOTSYS=/root-6.14.06/local -DARCH=haswell
