#!/bin/bash

sudo yum -y -t install epel-release
sudo yum repolist

sudo yum -y install cmake gcc gcc-c++ python python-devel python-pip geos geos-devel doxygen python-matplotlib python-sphinx

pip install --user setuptools
pip install --user cython
pip install --user numpy scipy shapely pint breathe sphinx-bootstrap-theme
pip install --user networkx nngt svg.path dxfgrabber

cd ..
mkdir build
cd build
cmake .. -Dwith-docs=ON

. ../set_dense_vars.sh

make && make install && make doc
