#!/bin/bash

sudo apt install cmake g++ python3 python3-dev python3-pip libgeos++-dev doxygen python3-matplotlib python3-tk libboost-dev

pip3 install --user setuptools
pip3 install --user "cython<3"
pip3 install --user numpy scipy pint pyneuroml
pip3 install --user --no-binary shapely, shapely
pip3 install --user sphinx breathe sphinx-bootstrap-theme
pip3 install --user networkx nngt svg.path dxfgrabber PyOpenGL

cd ..
mkdir build
cd build
cmake .. -Dwith-python=3 -Dwith-docs=OFF

. ../set_dense_vars.sh

make && make install
