To compile and run MinLinMo on Windows natively it's best to install Visual Studio with C++ support. It is now possible to install the vcpkg C++ library manager. To do this, open a Powershell and enter:

cd C:\
mkdir DEV
cd DEV

git clone https://github.com/microsoft/vcpkg.git

Go to the vcpkg folder by typeing:

cd vcpkg
.\bootstrap-vcpkg.bat
.\vcpkg integrate install

GSL can now be installed with vcpkg:
.\vcpkg install gsl gsl:x64-windows

Make a new C++ console project in Visual Studio and include the MinLinMo.cpp and MinLinMo.hpp files in the project and compile.

These instructions were first encountered here:
https://solarianprogrammer.com/2020/01/26/getting-started-gsl-gnu-scientific-library-windows-macos-linux/#gsl_installation_windows

After compilation, or to run the included executable file (GSL does not need to be installed), make sure gsl.dll and gslcblas.dll are in the same directory.

See also:https://zenodo.org/records/10149465
