#!/bin/bash

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

brew tap homebrew/core

PYVERSION="3.9"

echo 'export PATH="/Users/runner/Library/Python/$PYVERSION/bin:$PATH"' >> ~/.bash_profile

export PATH="/Users/runner/Library/Python/$PYVERSION/bin:$PATH"

brew install gcc@13 cmake "python@$PYVERSION" geos doxygen boost libomp

brew link gcc@13

export CXX=g++-13
export CC=gcc-13

pip3 install --user setuptools
pip3 install --user "cython<3"
pip3 install --user numpy scipy pint pyneuroml
pip3 install --user sphinx breathe sphinx-bootstrap-theme
pip3 install --user --no-binary shapely, shapely
pip3 install --user matplotlib networkx nngt svg.path dxfgrabber PyOpenGL

cd ..
mkdir build
cd build
cmake .. -Dwith-docs=ON -Dwith-python=3

make
make install
