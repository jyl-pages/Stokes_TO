# CMake生成脚本开发文档
在`simplex`项目中，我们使用3.18以上版本的CMake工具自动生成开发所需的Visual Studio工程。

3.18是因为CUDA_ARCHITECTURES功能需要3.18版本。

本文档介绍`simplex`使用的CMake生成脚本的基本情况。

`simplex`项目下一共有四个CMake脚本文件，都在`simplex/script`目录下，它们分别是：

* `functions.cmake`
* `common_func.cmake`
* `ext_config.cmake`
* `src_config.cmake`

下面我们具体介绍我们`CMake`脚本代码的组织方式。

# 生成脚本设计

`simplex`下一共有三种CMakeLists.txt文件，分别是：

* `simplex/proj/*`目录下的项目生成文件
* `simplex/src/*`目录下的`src`代码库生成文件
* 部分`simplex/ext/*`目录下的`ext`代码库生成文件

所有生成文件都会包含`functions.cmake`文件，并使用其中的接口。第一种生成文件会生成一个项目（project）和一个`.exe`生成目标（target）。后二者则只会生成一个`.lib`生成目标（target）。请注意辨析CMake下project和target的区别。一次CMake调用只会生成一个project，它对应于一个`Visual Studio`解决方案，但可以在这个project下生成若干个target，这些target可以是可执行文件，也可以是库文件（在`simplex`下，都是静态库）。

所有这些生成文件遵循一种统一的流程：

```cmake
init
use_ext_libs
use_src_libs
setup
```
具体参考[compilation-notes-zh.md](/compilation-notes-zh.md).

这个流程描述了我们设计生成脚本的基本思路。我们采用“先标记，后处理”的方案，对于每个target，记录五个LIST：

* 该target需要编译的所有源代码文件
* 该target需要用target_include_directories的所有目录
* 该target需要用target_link_libraries链接的所有二进制文件
* 该target依赖的所有target
* 该target需要被拷贝至可执行文件同一目录下的.dll文件

具体请参考`functions.cmake`里面的注释。它们均在`init_target()`中被创建。例如，第一个列表的名字叫做`${current_target}_srcs`，其中变量`current_target`储存当前target的名称。

为了自动处理依赖问题，我们在`add_subdirectory`后会自动propagate子项目的include目录和链接目录，具体见`common_func.cmake`下的`target_install_dep()`。这一点依赖于三组全局变量：`${current_target}_incs_global`、`${current_target}_libs_global`、`${current_target}_dlls_global`

# CMake脚本文件内容

`src_config.cmake`里面维护一个包含所有src代码库名称的LIST：`src_lib_lis`，还有`use_src_libs`，以及src代码库的CMakeLists.txt需要使用的setup_src_target宏。

`ext_config.cmake`和`src_config.cmake`差不多，但是它维护两个LIST，一个是所有只有头文件的第三方库列表`ext_header_lib_lis`，第二个是所有需要像`src`库那样编译成`.lib`的第三方库列表`ext_uncompiled_lib_lis`。

`functions.cmake`里面包含所有生成项目需要使用的公共函数。

`common_func.cmake`也包含一些公共函数，但是和`functions.cmake`相比更一般化。它和`functions.cmake`之间的界限略有些模糊，这是一个历史遗留问题。

# Linux编译的特殊处理

在Linux下，假如main.cpp引用了foo.a中的Foo()函数，而Foo()函数又引用了bar.a中的Bar()函数，则引用顺序应当为main.cpp、foo.a、bar.a。否则，因为编译器按顺序解析，如果bar.a在foo.a前面，则编译器只能从main.cpp中看到Foo()，从而不认为需要用Bar()，就会将其忽略，后面就会产生Bar()的undefined reference错误。因此，可以看到，common_func.cmake里面，merge_resources_from_target()中，如果系统是UNIX，就会把新的库加到前面而不是后面。在target_add_libs()中，也会在UNIX系统中，把新库加到前面。

# 开发人员列表

王梦迪
