# MinLinMo ARM-Linux (Raspberry pi)

To compile MinLinMo on an ARM processor with Raspbian or Ubuntu-type Linux, make sure the GSL library is installed:

_sudo apt install libgsl-dev_

Then compile with:

_g++ MinLinMo.cpp -lpthread -std=c++14  -O3 -lgsl -lblas -lpthread -march=armv8-a -o MinLinMo_
