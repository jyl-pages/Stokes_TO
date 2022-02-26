# 一、代码结构说明
simplex代码库分成三个部分：

* src代码库：位于simplex/src目录下的各文件夹，包括common、reservoir、geometry等。是由实验室维护的，提供编写代码使用的基础数据结构和算法等。
* ext代码库：位于simplex/ext目录下的各文件夹，包括eigen、nanoflann等。是项目代码或src代码需要用到的第三方库。
* proj工程：位于simplex/proj目录下的各文件夹，包括fluid_euler、levelset等。是每一个具体项目的代码。

编译的时候，以proj代码目录为基础，它会根据设置，自动包含或编译相应的src/proj代码库。

# 二、proj工程CMakeLists.txt编写说明

为了方便，我们在下文中用`my_proj`指代工程名称。它的`CMakeLists.txt`的格式如下：

```cmake
cmake_minimum_required(VERSION 3.18)
set(proj_name my_proj)
project(${proj_name})
include(${PROJECT_SOURCE_DIR}/../../../simplex/script/functions.cmake)
init_project(${PROJECT_SOURCE_DIR}/../../../simplex)
use_ext_libs(LINK <list-of-ext-libs>)
use_src_libs(SELF <list-of-src-libs>)
setup_exe_target()
```

你需要更改的只有三条命令：

* `set(proj_name my_proj)`
* `use_ext_libs(LINK <list-of-ext-libs>)`
* `use_src_libs(SELF <list-of-src-libs>)`

第一条只需要把`my_proj`替换成你的项目名称即可。后面的两条，对于大部分模拟器而言，如下设置就够用了：

```cmake
use_ext_libs(LINK "cuda;openmp;physbam;")
use_src_libs(SELF "physics;")
```

当然，你可以自定义修改它们。例如，如果你希望在Visual Studio工程文件中增加某个宏定义，则应当在`setup_exe_target()`之前使用`add_definitions`指令，如：

`add_definitions(-DUSE_SOME_LIB)`

下面我们具体介绍`use_ext_libs`和`use_src_libs`这两条指令的用法。

### 使用ext代码库
这条指令设置`my_proj`使用的第三方代码库，例：
```cmake
use_ext_libs(LINK "cuda;openmp;physbam;")
```
所有ext库均使用直接链接`.lib`文件的模式，`my_proj`永远不会尝试编译它们。该指令的第一个参数固定为`LINK`，意义就是采用这种链接编译模式。第二个参数是一个`CMake`语言中的list，它是一个由分号隔开的字符串。包含所有你想在本项目中使用的所有第三方库。在这个例子中，就是`physbam`、`openmp`和`cuda`.

注意，（在src代码库的`SELF`也就是默认模式下）项目会继承你在这些src代码库中使用的所有第三方库。所以你不需要考虑src中都用了什么第三方库。例如，只要你使用了`src/common`，你就永远不需要在这里添加`eigen`一项，因为`common`包含了`eigen`.

你可以参考`simplex/ext`目录下查看有哪些可以使用的第三方库。注意还有`openmp`和`cuda`。在这里使用`cuda`意味着打开CUDA语言支持，不使用则关闭支持。你不需要做其他任何操作，生成脚本会自动处理一切事情。

### 使用src代码库
这条指令设置`my_proj`使用的src代码库，你可以在`simplex/src`中找到它们，它们是一些你编写模拟程序所需要的基础数据结构和算法。例：

```cmake
use_src_libs(SELF "viewer;physics")
```

第一个参数指明使用src代码库的模式。<font color='red'> 目前只全面支持`SELF`模式 </font>，还有一个相关不完整的`COMMON`模式。

* 如果设置为`SELF`，则`my_proj`工程将在`simplex/bin/my_proj`目录下生成所有用到的src代码库对应的.lib文件，并链接它们以生成可执行文件。不同的工程项目之间的编译过程不会相互干扰。这样，如果你修改了`simplex/src`下（且被实例化至同目录下的.cpp文件）的程序代码，`my_proj`项目会自动检测到这一变化，并重新编译相应的.lib文件。
* 如果设置为`COMMON`，则`my_proj`工程将不会试图编译src代码库，而只是链接已编译好的，位于`simplex/lib/simplex`目录下的src代码库对应的.lib文件。这些文件需要通过编译`simplex/proj/_install`工程生成。这样，如果你修改了`simplex/src`下的程序代码，必须重新编译`simplex/proj/_install`，才能在`my_proj`中使用新的代码。

第二个参数，和`use_ext_libs`类似，包含你所需要使用的一切`src`库名称。正如刚才所说，对于大部分模拟器，你只需要将其设置为`physics`即可。

在推荐的`SELF`模式下，我们的编译指令可以自动解决一切互相依赖问题。例如，库`reservoir`依赖于库`common`，所以如果你在这里使用了`reservoir`，代码就会自动检测到这一点，自动编译、链接`common`，并自动添加`common`所需的的一切支持（例如，使用`eigen`库）。你可以看到，在上面的示例中，我们只使用了`viewer`和`physics`两个库，在这种情况下，生成脚本会自动检测依赖关系，并把`common`、`reservoir`、`geometry`、`solver`……等加入代码。

特别地，你不必担心关于CUDA的任何问题。如果你在`use_ext_libs`中设置了`cuda`，`physics`库将自动添加`parallel`库。如果不设置`cuda`，生成脚本将检测到这一点，略去`parallel`库的生成。

当然，如果你没有使用`physics`而又想要使用CUDA，则可以在`use_ext_libs`中设置`cuda`，并且在这里添加`parallel`。如果你又想编译一个纯CPU版本，则只需要简单地在`use_ext_libs`的参数中删去`cuda`，生成脚本将自动忽略`parallel`，不用手动删除。

但是`COMMON`模式尚不支持这些功能。<font color='red'> TODO:支持`COMMON`模式 </font>

## 屏蔽代码模块

目前，我们的代码有两套GPU求解器，第一套是`src/parallel`下的传统求解器（下称parallel solver），第二套是`src/common_cpx`和`src/poisson_cpx`下的线性系统求解器（下称cpx solver）。目前，如果你在`ext`中使用了`cuda`库，同时又在`src`中使用了`physics`库，那么代码会同时编译parallel和cpx两套功能相同的求解器，这可能是冗余的。我们提供了`ignore_src_modules`函数用来关闭求解器的编译。你只需要在`use_src_libs`语句前使用它：

```cmake
ignore_src_modules("cpx;")
use_src_libs(...)
```

这样就能让代码不编译cpx solver。换成`parallel`就不编译parallel solver，换成`parallel;cpx`就两个都不编译。

# 三、向simplex中新增一个src库

为了方便起见，我们在下文中用`my_src`指代这个src代码库的名称。你需要做三件事情。

1. 在`simplex/src/my_src`目录下添加所有代码文件。
1. 在`simplex/script/src_config.cmake`的第一行`set`语句中添加：`set(src_lib_lis ... my_src)`。
1. 编写`simplex/src/my_src/CMakeLists.txt`文件

### src代码库CMakeLists.txt格式
`simplex/src/my_src/CMakeLists.txt`文件遵循固定格式如下：
```cmake
init_target(<name-of-src-lib>)
use_ext_libs(COMPILE <list-of-ext-libs>)
use_src_libs(COMPILE <list-of-src-libs>)
setup_src_target()
```
第一行的含义是设置这个代码库的名称，和它所在的文件夹相同，如`common`、`geometry`。对于`my_src`就是：

    init_target(my_src)
第二行的含义是设置所有需要使用的ext代码库。例如：

    use_ext_libs(COMPILE "ipopt;nlopt;mma")
它和`proj`工程CMakeLists.txt里面的`use_ext_libs`基本相同。唯一的不同是，第一个参数也是`COMPILE`而非`LINK`。

第三行设置`my_src`依赖的所有src代码库。注意也是COMPILE模式。和`proj`相同，你也只需要添加直接依赖的库。例如：

    use_src_libs(COMPILE "solver;geometry;parallel;")

第四行是固定的，无需更改。

# 四、向simplex中新增一个第三方库

假设你在项目中使用了一个名叫`my_ext`的第三方库，并且希望将它添加到simplex代码库中。在这种情况下，你需要做两件事情：

1. 在`simplex/ext`目录下添加相应的目录`simplex/ext/my_ext`，并向其中添加这个库包含的所有文件。
2. 修改`simplex/script/ext_config.cmake`文件，以支持编译`my_ext`库。

具体的修改方式，根据`my_ext`库的特点，分成三种情况：

* `my_ext`中仅包含头文件(.h和.hpp)
* `my_ext`由一系列头文件(.h,.hpp)和源代码文件(.cpp)组成，我们需要手动将其编译成库文件。
* 其他情况。

### 情况1：`my_ext`仅包含头文件

在这种情况下，只需要在项目中包含`my_ext`的所有文件，就能正确地使用这个库。此时，请在`ext_config.cmake`文件开头的

    set(ext_header_lib_lis ...)
命令中添加`my_ext`项。例如：

    set(ext_header_lib_lis amgcl ... perlin_noise my_ext)
这个命令的含义是创建一个名叫`ext_header_lib_lis`的`LIST`对象。随后，CMake代码就会自动判断何时将`my_lib`加入项目。

### 情况2：`my_ext`需要手动编译
如果`my_ext`的目录结构只是在`simplex/ext/my_ext`下有一些`.h,.hpp,.cpp`文件，不含任何子目录，（也就是说，和`simplex/src`下的各个目录结构相同），那么我们需要运行`simplex/proj/_install`项目，将它编译成静态库，把相应的`.lib`（或UNIX下的`.a`）文件放在`simplex/lib/ext`下。随后，在项目中使用时，`use_ext_libs()`指令会自动链接这些文件。

此时，你需要在`simplex/ext/my_ext`目录下添加一个`CMakeLists.txt`文件。它的内容是：

```cmake
init_target(<name-of-src-lib>)
use_ext_libs(COMPILE <list-of-ext-libs>)
use_src_libs(COMPILE <list-of-src-libs>)
setup_exe_target()
```

这和`src`库的CMakeLists.txt格式完全相同，只是最后一行是`setup_ext_target`而非`setup_src_target`。以`tiny_obj_loader`为例：

```cmake
init_target(tiny_obj_loader)
#use_ext_libs(COMPILE)
use_src_libs(COMPILE "reservoir;")
setup_ext_target()
```
随后，和情况1类似，你只需要在`ext_config.cmake`文件开头的

    set(ext_uncompiled_lib_lis ...)
命令中添加`my_ext`即可。注意需要修改的命令不是情况1的那条命令。

### 情况3：其他情况

如果安装`my_ext`的过程比较复杂，不在前两种情况之内，你就需要更改`ext_config.cmake`文件中的`use_ext_libs()`宏。你可以参考其中已有的代码，例如支持使用`libtorch`库的代码：

```cmake
if("libtorch" IN_LIST lib_lis)
  message(STATUS "Use ext lib: libtorch")
  set(CMAKE_PREFIX_PATH ${simplex_path}/ext/libtorch/share/cmake/Torch)
  find_package(Torch REQUIRED)
  list(APPEND ${current_target}_libs ${TORCH_LIBRARIES})
endif()
```

这里可能涉及到一些`CMake`生成脚本的内部结构。比如`${current_target}_libs`是一个`LIST`类型的变量，指代当前编译目标所需要链接的所有库文件。你可能需要阅读`simplex/script`下的`functions.cmake`、`ext_config.cmake`、`src_config.cmake`、`common_func.cmake`文件。如果遇到任何问题，或感到代码难以理解，请直接询问`simplex`代码库维护人员。

# 作者列表

王梦迪
