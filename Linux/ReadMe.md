If the GNU Scientific Library is not installed you will need the following libraries:

libgsl.so.27
libgslcblas.so.0 (could be a soft link to libgslcblas.so.0.0.0)

On Ubuntu/MINT type Linux they should be put in the following directory:
/usr/lib/x86_64-linux-gnu/

Compile MinLinMo for Linux with:
g++ MinLinMo.cpp -march=core-avx2 -lgsl -lblas -lpthread -O3 -o MinLinMo

The GNU Scientific library can easily be installed on Ubuntu/MINT:
sudo apt-get install libgsl-dev
