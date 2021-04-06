cd ..
mkdir build
cd build

set DISTUTILS_USE_SDK 1

cmake .. -DCMAKE_GENERATOR_PLATFORM=x64
cmake --build . --config Release --target INSTALL
