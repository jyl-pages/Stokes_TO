## Linking OpenCV Libraries to Simplex

-- This document was created by Yijun Li.

### Prerequisites

1. GCC
2. CMake
3. OpenCV

* For Linux users, a higher version of GCC is required to support C++17 standard.
* If root operations are allowed, try `apt-get` first before building OpenCV on the machine on from source. Please refer to official document for more help. This doc will focus on setup w/o root.
* When using CMake to build OpenCV, probably 3.5 is required. Add the higher CMake to your path or manually call it in terminal.
* Remember to configure `CMAKE_INSTALL_PREFIX` because this will be the ultimate path CMake is looking for in your own project. Here I will set it to `/path-to-opencv/opencv-4.2.0/build/install`. Without modifying this parameter, the old version will be found and may cause errors like `undefined variable`.
* As explained in the [OpenCV Official Guide](https://docs.opencv.org/4.2.0/d7/d9f/tutorial_linux_install.html), make sure to `mkdir` and `cd` the `/build` directory, and configure a different directory like `/build/install` in step 2. No need to `mkdir install`. Run `make -j7 ` and `make install` in the `build` dir.



### Linking OpenCV to your project

1. Compile `/_install` as usual
2. Modify the `CMakeList.txt`:
   
    First add after the cmake min version, before your `project(${proj_name})`:
    >SET(OpenCV_DIR "/path-to-opencv/opencv-4.2.0/build/install")   
    SET(CMAKE_C_COMPILER "/data00/home/liyijun.amber/gcc-8.2.0/objdir/bin/gcc")
    SET(CMAKE_CXX_COMPILER "/data00/home/liyijun.amber/gcc-8.2.0/objdir/bin/g++")
    find_package( OpenCV REQUIRED )

    Then add after `link_project()`:
    > target_link_libraries( ${proj_name} ${OpenCV_LIBS} )

    These messages can also be added in the end to debug:
    > message(STATUS "OpenCV library status:")
    message(STATUS "    config: ${OpenCV_DIR}")
    message(STATUS "    version: ${OpenCV_VERSION}")
    message(STATUS "    link libraries: ${OpenCV_LIBS}")
    message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")


### Reference:
[OpenCV Official Guide](https://docs.opencv.org/4.2.0/d7/d9f/tutorial_linux_install.html)

[Guid for nm command to check symbols in Linux](https://www.ibm.com/support/knowledgecenter/ssw_aix_71/n_commands/nm.html)