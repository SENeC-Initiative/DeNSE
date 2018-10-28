mkdir build
cd build

CC=gcc CXX=g++ cmake .. -Dwith-docs=ON

. ../set_dense_vars.sh

make
CC=gcc CXX=g++ make install
