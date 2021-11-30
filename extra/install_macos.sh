#!/bin/bash

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

brew tap homebrew/core

brew install gcc@10 cmake python3 geos doxygen boost libomp

brew link gcc@10

pip3 install setuptools
pip3 install cython
pip3 install numpy scipy pint pyneuroml
pip3 install sphinx breathe sphinx-bootstrap-theme
pip3 install --no-binary shapely, shapely
pip3 install matplotlib networkx nngt svg.path dxfgrabber PyOpenGL

cd ..
mkdir build
cd build
CC=gcc-10 CXX=g++-10 cmake .. -Dwith-docs=ON -Dwith-python=3

make
CC=gcc-10 CXX=g++-10 make install
