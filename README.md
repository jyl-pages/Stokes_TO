
GPU_Stokes_TO
================
Dependencies
--------------------
CUDA 11.0+

Compilation
-----------------
Compilation has been tested with Microsoft Visual Studio 19 (under Windows 10).

To compile libWetCloth, you'll need CMake-GUI 3.21+ (https://cmake.org) on Windows.

Step 1. Compile dependencies: In CMake-GUI, browse source code to .../simplex/proj/_install; browse build to .../simplex/build/_install/cmake-build-win; configue with VS16 2019 and default compiler; generate and open project. Turn to Release mode and build, then close the project.

Step 2. Compile opengl viewer: Repeat Step 1 by replacing "_install" with "opengl_viewer".

Step 3. Compile project: In CMake_GUI, browse source to .../complex/proj/PainlessSolver/proj/fluid_topo; browse build to .../complex/build/fluid_topo/cmake-build-win; configure, generate and open project. Turn to Release mode and build.

Step 4. Copy files under .../simplex/bin/win/opengl_viewer/Release to .../complex/bin/win/fluid_topo/Release.

Run the Demo
--------------------
In the command tool, go to the .../complex/bin/win/fluid_topo/Release folder, type:
```
./fluid_topo.exe -driver 0 -test 0 -s 256 -o output -frac 0.15 -lf 40
```
Here, -s for resolution, -o for file name to write result, -frac for volume fraction, and -lf for optimization it

