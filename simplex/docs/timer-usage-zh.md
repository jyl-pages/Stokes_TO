# 使用Timer类测量代码运行时间
模拟代码通常对性能有较高的要求。common/Timer.h中的Timer类提供了记录代码运行时间的功能。

它的用法很简单：假设你的主循环由三个部分A,B,C组成。你只需要创建一个Timer类（不妨叫它timer），然后在循环的开头调用
```c++
timer.Begin_Loop();
```
在A,B,C三部分结束时依次调用

```c++
timer.Record("A");
timer.Record("B");
timer.Record("C");
```

然后在循环结束时调用

```c++
timer.End_Loop_And_Output(std::cout);
```

代码就会依次输出A,B,C三部分的本次循环用时和平均用时（单位是秒）。最后的std::cout也可以换成一个ofstream，这样就会把运行时间记录在文件中。

# 作者列表

王梦迪
