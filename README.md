
GPU_Stokes_TO
================
Dependencies
--------------------
CUDA 11.0+

Compilation
-----------------
Compilation has been tested with Microsoft Visual Studio 19 (under Windows 10).

To compile libWetCloth, you'll need CMake-GUI 3.21+ (https://cmake.org) on Windows.



Run the Demo
--------------------
To run the demo of libWetCloth, you may simply use the command line argument *-s [scene_file]* to specify the scene to be loaded. Besides of the WetCloth library, we also have a user interface named WetClothApp where you may watch and tune the simulation interactively. For example, you may type
```
./WetClothApp -s ../assets/unit_tests/simple_cloth.xml
```
under the `<build>/libWetCloth/App` directory to run the simulation of the scene containing a water ball splashes on a small cloth. 

All the parameters can be modified offline in the scene description XML files. Some can be changed online in the user interface provided by the demo program.

USAGE: 
```
   ./WetClothApp -s <string> [-i <string>] [-o <integer>] [-g <integer>] [-d <boolean>] [-p <boolean>] [--] [--version] [-h]

Where: 
   -s <string>,  --scene <string>
     (required)  Simulation to run; an xml scene file

   -i <string>,  --inputfile <string>
     Binary file to load simulation pos from

   -o <integer>,  --outputfile <integer>
     Binary file to save simulation state to

   -g <integer>,  --generate <integer>
     Generate PNG if 1, not if 0

   -d <boolean>,  --display <boolean>
     Run the simulation with display enabled if 1, without if 0

   -p <boolean>,  --paused <boolean>
     Begin the simulation paused if 1, running if 0

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

Surface Reconstruction and Rendering with Houdini
--------------------------------------------------------
The Houdini projects are also provided in the "houdini" folder, which are used for surface reconstruction and rendering purposes. Our simulator can generate data that can be read back by the Python script in our Houdini projects.

You may run the demo with the "-o" option, *under the folder containing Houdini projects* (by default, [project source]/houdini/). For example, you may type
```
../build/libWetCloth/App/WetClothApp -s ../assets/parameter_tests/pore_test_mid.xml -o 20
```
under the folder *containing Houdini projects* to generate data per 20 time steps (for this example we use 0.0002s for the simulation time step, and 0.004s for the rendering time step. Hence 20 time steps is taken for the data generation). The simulation code will create a folder with the name of the scene file under this folder ([project source]/houdini/pore_test_mid for this example), which can be read back by the Houdini project with the same name.

After some data generated, you may open the corresponding Houdini project to watch and operate on them. We use the nodes with suffix "bake" to indicate the usage of baking. For example in the "fluids_bake" node, the fluid particles will be read and used to reconstruct a polygonal liquid surface, which will then be stored as Houdini geometry files.

The baked geometries are then read back by the nodes with prefix "geo", which can be used for rendering. 

You may check the Mantra nodes in the output classifier for the rendering settings. You can also connect these nodes to HQueue nodes for distributed rendering on a render farm.

Houdini 16.5 is used for our case. You may need an equivalent or higher version to open the project files. For more information and tutorials of Houdini, please visit the SideFX website (https://www.sidefx.com/).

Contact
-----------
Please contact the author (yf2320@columbia.edu) for questions and bug report, or consulting for the usage of this code base.

BibTex Citation
----------------------
@article{Fei:2018:MMS:3197517.3201392  
 author = {Fei, Yun (Raymond) and Batty, Christopher and Grinspun, Eitan and Zheng, Changxi},  
 title = {A Multi-scale Model for Simulating Liquid-fabric Interactions},  
 journal = {ACM Trans. Graph.},  
 issue_date = {Aug 2018},  
 volume = {37},  
 number = {4},  
 month = aug,  
 year = {2018},  
 pages = {51:1--51:16},  
 articleno = {51},  
 numpages = {16},  
 url = {http://doi.acm.org/10.1145/3197517.3201392},  
 doi = {10.1145/3197517.3201392},  
 acmid = {3201392},  
 publisher = {ACM},  
 address = {New York, NY, USA},  
 keywords = {fluid dynamics, particle-in-cell, mixture theory, two-way coupling, wet cloth},  
}  
