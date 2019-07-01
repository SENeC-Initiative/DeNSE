#!/bin/bash

sudo apt install cmake g++ python python-dev python-pip libgeos++-dev doxygen python-matplotlib python-tk libboost-dev

pip install --user setuptools
pip install --user cython
pip install --user numpy scipy shapely pint sphinx breathe sphinx-bootstrap-theme
pip install --user networkx nngt svg.path dxfgrabber PyOpenGL

cd ..
mkdir build
cd build
cmake .. -Dwith-docs=ON

. ../set_dense_vars.sh

make && make install && make doc
