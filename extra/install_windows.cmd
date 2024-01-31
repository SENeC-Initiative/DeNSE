IF "%1"=="" ( SET "pyversion=3" ) ELSE ( SET "pyversion=%1" )

cd ..
mkdir build
cd build

cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -Dwith-python=%pyversion%
cmake --build . --config Release --target INSTALL
