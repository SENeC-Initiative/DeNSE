#!/bin/bash

sudo pacman -S cmake gcc python python-pip geos boost boost-libs doxygen python-matplotlib python-sphinx

pip install --user setuptools
pip install --user cython
pip install --user numpy scipy shapely pint breathe sphinx-bootstrap-theme
pip install --user networkx nngt svg.path dxfgrabber PyOpenGL

mkdir build
cd build
cmake .. -Dwith-docs=ON

. ../set_dense_vars.sh

make && make install && make doc
