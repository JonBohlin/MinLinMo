If the GNU Scientific Library is not installed you will need the following libraries:

libgsl.so.27 libgslcblas.so.0 (could be a soft link to libgslcblas.so.0.0.0)

Compile MinLinMo for OS X with Apple silicon using (-I and -L point to GSL include and library files respectively):

g++ MinLinMo.cpp -lpthread -std=c++17 -o MinLinMo -O3 -lgsl -lblas -I /opt/homebrew/Cellar/gsl/2.7.1/include/ -L /opt/homebrew/Cellar/gsl/2.7.1/lib/ -march=armv8-a

The GNU Scientific library can easily be installed on OS X with Homebrew: brew install gsl
