cd ..
mkdir build
cd build

set DISTUTILS_USE_SDK 1

cmake ..
cmake --build . --config Release --target INSTALL
