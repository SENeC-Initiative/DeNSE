#!/bin/bash

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

brew tap homebrew/core

brew install gcc-8 cmake python3 geos doxygen boost

pip3 install setuptools
pip3 install cython
pip3 install numpy scipy pint breathe sphinx-bootstrap-theme
pip3 install --no-binary shapely, shapely
pip3 install matplotlib networkx nngt svg.path dxfgrabber PyOpenGL

cd ..
mkdir build
cd build
CC=gcc-8 CXX=g++-8 cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local -Dwith-docs=ON -Dwith-python=3

make
CC=gcc-8 CXX=g++-8 make install
make doc
