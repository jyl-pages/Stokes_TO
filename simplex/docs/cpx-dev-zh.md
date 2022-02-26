# CPX求解器开发文档

`simplex/src`目录下的`common_cpx`和`poisson_cpx`两个目录为我们提供了一个可以求解$512^3$级别计算域的，基于GPU的高性能泊松方程组求解器。我们将其称为CPX求解器。接下来，我们介绍它的基本情况。

# 共轭梯度法

`common_cpx/ConjugateGradient.h`提供了预优共轭梯度（PCG）的函数。

PCG是新引入一个对称正定矩阵$\bm{C}$，然后把求解方程$\bm{Ax}=\bm{b}$的问题变成

$$\bm{\tilde{A}}\bm{\tilde{x}}=\bm{\tilde{b}}.$$

这个$\bm{\tilde{A}}=\bm{C}^{-1}\bm{A}\bm{C}^{-1}$，$\bm{\tilde{x}}=\bm{Cx}$，$\bm{\tilde{b}}=\bm{C}^{-1}\bm{b}$.

然后在求解中，对CG改动一下，在第i步计算$\bm{z}_i=\bm{M}^{-1}\bm{r}_i$，这个$\bm{r}_i=\bm{b}-\bm{A}\bm{x}_i$是第$i$步残差，其中$\bm{x}_i$是第$i$步得到的结果。


# 开发人员列表

冼臧越洋、刘瑾源、王梦迪
