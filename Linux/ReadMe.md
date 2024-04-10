# Compiling MinLinMo on Linux

MinLinMo depends on the GNU Scientific library (GSL). GSL can easily be installed on an Ubuntu-type Linux:

_sudo apt-get install libgsl-dev_

If the GNU Scientific Library is not installed you will need the following libraries:

_libgsl.so.27_

_libgslcblas.so.0_ (can be a soft link to libgslcblas.so.0.0.0)

On an Ubuntu type Linux the library files should be put in the following directory:

_/usr/lib/x86_64-linux-gnu/_

MinLinMo can now be compiled for Linux with:

_g++ MinLinMo.cpp -march=core-avx2 -lgsl -lblas -lpthread -O3 -o MinLinMo_

**Compiling MinLinMo for Linux on ARM based CPUs (Raspberry Pi, for instance) see separate ARM directory.**
