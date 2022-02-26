# Use Vcpkg to Install Libraries

You can use Vcpkg to install and manage third-party C++ libraries. It can free you from the burden of libraries in simplex\ext.

We assume that you use Windows system in this tutorial, for other operator systems and full documentation, see:
https://github.com/microsoft/vcpkg

### Install Vcpkg

1. Assign an install path for Vcpkg. Here we take C:\Program Files :

        $ cd "C:\Program Files"
2. Download Vcpkg from GitHub:

        $ git clone https://github.com/microsoft/vcpkg.git
3. Install Vcpkg with script:

        $ .\vcpkg\bootstrap-vcpkg.bat

### Install C++ packages with vcpkg
Take Eigen as example here.

1. Navigate to Vcpkg path:

        $ cd "C:\Program Files\vcpkg"
2. Search the package which you want to install (Eigen here)

        $ .\vcpkg.exe search Eigen
3. You will see a list of packages like:

        ceres[eigensparse]                    Use of Eigen as a sparse linear algebra library in Ceres
        eigen3               3.3.7-5          C++ template library for linear algebra: matrices, vectors, numerical solvers,...
        libigl               2.1.0-2          libigl is a simple C++ geometry processing library. We have a wide functionali...

    Find the package name in this list (eigen3 here)
4. Install third-party package:

        $ .\vcpkg.exe install eigen3:x64-windows
    "x64-windows" here represents your OS. If you want to compile 32-bit version, then it's "x86-windows".
5. Check if the package is correctly installed:

        $ .\vcpkg.exe list
    In case of success, you will find womething like:

        eigen3:x64-windows                                 3.3.7-5          C++ template library for linear algebra: matrice...
6. Integrate the newly-installed package into Visual Studio:

        $ .\vcpkg.exe integrate install
7. Open CMake GUI, "Configure" and "Generate" your project again, like previously mentioned in this tutorial.
    Now open your project in Visual Studio. You can use Eigen library now!
