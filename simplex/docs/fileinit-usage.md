# 设置模拟器代码的参数
一份模拟代码往往需要设置一些常数，例如：Eulerian网格的分辨率和尺寸、所模拟流体的密度、粘性系数、固体的杨氏模量等等。有时候，你需要花费大量时间调整这些参数，以保证代码达到你想要的效果。

设置这些参数通常有两种办法。

第一种，是在相应`Driver`类的`Initialize()`函数中设置。这样做的好处是：
* 可以对多个实验配置不同的参数，例如你或许想在实验1中模拟粘度为`0.001`的水，而在实验2中模拟粘度为`10`的蜂蜜。
* 参数设置是代码的一部分，修改记录会被`git`自动保存。

但是这种办法有一个缺点，即它完全是在代码内部的，每次修改都必须重新编译，耗费大量编译时间，也不利于同时运行多组参数以搜索最佳方案。

当然，我们可以不需要编译而修改参数，也就是第二种：

在`main`函数中使用`ParseArgs`类，通过命令行传入相应参数。

这样做的好处是可以在不重新编译的情况下，只通过修改命令行参数的方式修改代码的运行模式，但也有如下缺点：

* 这些运行参数并不在代码中出现，必须手工记录，否则会丢失。
* 一些模拟器拥有数量较多的参数，这会让命令行变得极为冗长，而且容易遗漏或错误设置变量。
* 每次添加或删除常数都必须修改`main.cpp`文件。
* 只能设置一组默认值，不方便有多组实验的情况。

为了解决这些缺点，我们在`src/common/FileInit.h`中提供了一种从文件中读取参数的方法，它还包含强大的自动查错功能，接下来我们描述它的用法。

# FileInit的用法

## 1. 项目代码中的准备工作
为了使用`FileInit.h`中的接口，你需要在`SimulatorParams`的继承类中存放参数。例如，假设你需要写一个流体-固体耦合模拟器，其中包含3个参数：

* 流体的粘性
* 流体的密度
* 固体的弹性

你可以用两个`SimulatorParams`的继承类存放这三个参数：

```c++
class FluidParams: public SimulatorParams{
public:
  real viscosity;
  real density;
};
class SolidParams: public SimulatorParams{
public:
  real elasticity;
};
```

你需要在每个类中做两件事情：
* 给它起一个独特的标识符，并保证基类中的变量`std::string name`永远被置为这个标识符。设置的方法可以是通过基类`SimulatorParams`的构造函数`SimulatorParams(std::string)`。标识符中不得含有空格或字符`#@`，这两个字符是保留字符。例如，你可以把`FluidParams`的标识符记为`fluid`，把`SolidParams`的标识符记为`solid`。如果你忘了设置，它将是`SimulatorParams`类中赋予的默认值`""`（空字符串），并在后面读取的时候报错，所以无需过度担心。
* 实现`SimulatorParams`类的虚函数`virtual void Read(std::istream& input)`.这个函数的作用是从一个`std::istream`读取类中的所有参数。由于有时候你可能只想写一个参数类而不想使用文件读取功能，因此我们没有把`SimulatorParams`中的该函数设为空虚函数（从而让编译器逼迫你实现它），但是如果你没有在继承类中实现它，而又试图读取这个继承类，程序就会自动执行`SimulatorParams`中定义的该函数，它会输出错误信息，提示你实现这个类，并直接退出。

示例代码：

```c++
class FluidParams: public SimulatorParams{
public:
  FluidParams():SimulatorParams("fluid"){}
  real viscosity;
  real density;
  virtual void Read(std::istream& input) {
    input >> viscosity >> density;
  }
};
class SolidParams: public SimulatorParams{
public:
  SolidParams():SimulatorParams("solid"){}
  real elasticity;
  virtual void Read(std::istream& input) {
    input >> elasticity;
  }
};
```

## 2.参数文件的格式

接下来，你需要编写一个参数文件。它的格式如下：

```
# Some descriptions
@<param_name1>
[value1_1]

[value1_2] # some comments
...
@<param_name2>
[value2_1] [value2_2]
[value2_3]
...
```

也就是说：
* 每一个参数块以`@+参数类标识符`为开头，紧跟一行或多行内容。这些内容应当能被相应参数类的`Read`函数读取。
* 你可以像在`Python`中一样，在任意位置使用单行注释符号`#`。我们的代码会自动忽略注释。
* 在参数文件中不影响语义的地方，你可以任意添加空行或空格。例如，例子中`[value1_1]`和`[value1_2]`中间的空行，它会被我们的代码自动忽略。
* 参数块的顺序不会产生影响，你可以任意布置。但是各参数块的标识符不能重复，也不能为空。否则我们的代码会检测到这一点，并报错退出。唯一的例外是，后面不跟随任何内容的`@`会被忽略，它不包含任何信息。

例如：
```
#Parameters of my fluid-solid coupling simulator
@fluid
1e3 # density of FluidParams
0.01#viscosity

@solid
#just elasticity
100.
```

## 3. 读取参数文件
使用`FileInit::Read_Params_From_File()`函数，从文件中读取各个参数类。例如，假设参数文件的路径为`params.txt`，则代码为：

```c++
FluidParams fluid_params;
SolidParams solid_params;
FileInit::Read_Params_From_File("params.txt",&fluid_params,&solid_params);
```

注意事项如下：
* `FileInit::Read_Params_From_File()`函数使用了可变参数列表技术，就像`printf()`一样，你可以传入任意多个`SimulatorParams`继承类的指针。
* 该函数基于`std::map`查询对应的参数块，因此参数文件中各块的顺序，和函数调用列表中的顺序均不影响，可以任意布置。
* 我们的代码会检测每个参数块是否被没有错误地读完，否则将报错退出，并输出相应的调试信息。例如，参数个数过少、过多，或本应是数字的地方有字母插入，代码都会报错。

# 作者列表

王梦迪
