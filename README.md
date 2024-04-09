# MinLinMo
Fast and efficient (n<<p) variable selection for linear models

## Compiling and installing MinLinMo
MinLinMo was written in C++ version 14. It has been developed for both Intel and ARM
processors. For Intel processors both multi-threading and AVX2 are employed for maximum
performance. For ARM processors, NEON is used instead. As AVX2 supports 256 bit registers parts
of MinLinMo will likely run faster on Intel processors. During the model building phase, which is
also the most time consuming phase, multi-threading was not found to improve performance.

To compile MinLinMo for Ubuntu-type Linux make sure that the GNU Scientific Library (GSL) is
installed. This can be installed with:

 _sudo apt install libgsl-dev_

MinLinMo can now be compiled with:

 _g++ MinLinMo.cpp -march=core-avx2 -lgsl -lblas -lpthread -O3 -o MinLinMo -std=c++14_ 

For Apple silicon, use Homebrew to install GSL:

 _brew install gsl_

To compile MinLinMo for the ARM (Apple silicon) OS X version use:

 _g++ MinLinMo.cpp -lpthread -std=c++14 -O3 -lgsl -lblas -I /opt/homebrew/Cellar/gsl/2.7.1/include/ -L /opt/homebrew/Cellar/gsl/2.7.1/lib/ -march=armv8-a -o MinLinMo_

For ARM Linux the following should hopefully work, providing the gsl is installed (but has not been tested):

 _g++ MinLinMo.cpp -march=armv8-a -lgsl -lblas -lpthread -O3 -o MinLinMo -std=c++14_

For Intel OS X this will hopefully work:

 _g++ MinLinMo.cpp -lpthread -std=c++14 -O3 -lgsl -lblas -I /opt/homebrew/Cellar/gsl/2.7.1/include/ -L /opt/homebrew/Cellar/gsl/2.7.1/lib/ -march= core-avx2_

For Windows (and Intel-based Macs with Boot camp), the easiest way to run MinLinMo is to downloaded the .exe and .dll files and make sure the files are in the same directory.

To compile MinLinMo from scratch, it will likely be easiest to install GSL through Visual Studio 2017 or later. This can be done
as follows:

- Create a new C++ project, right click the project name under Solution Explorer and choose
“Manage Nuget Packages”
- Click the Browser tab and search with the keywords “Microsoft.gsl”, which will filter out
the Microsoft GSL version. It can now be installed by clicking the install button

MinLinMo can now be compiled with the Microsoft C++ compiler from Visual Studio.

To test MinLinMo, first download the included _mtcars_ dataset from the R packacge, i.e. the files mpg.txt (outcome) and mtcars.csv (predictors). For Linux and OS X, run MinLinMo from a path to the directory it was installed or in the current directory:

 _./MinLinMo mpg.txt mtcars.csv_

 For Windows, the dll files must be in the same directory as the MinLinMo.exe file. MinLinMo can then be run with a path to the directory it was installed, or in the current by typing:

 _MinLinMo.exe mpg.txt mtcars.csv_
 

 

