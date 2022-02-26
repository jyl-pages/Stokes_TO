All our kernels are isotropic ands follow the form

$$
W_d(\vec{r})=\frac{\alpha_d}{h^d}\phi\left(\frac{|\vec{r}|}{h}\right).
$$

Where $d\in \{1,2,3\}$ is the dimension number. $\phi(x)$ is some scalar weight function truncated at $1$, which means $\forall x\geq 1, \phi(x)=0$, so $W_d(\vec{r})$ is truncated at a variable radius $h$. And $\alpha_d$ is a constant coefficient related with $\phi$ and dimension number $d$. You can verify that for all dimensions, when $\int\alpha_d\phi=1$ is satisfied, $\int W_d(\vec{r})=1$ still holds.

And we can calculate the gradient of $W_d$ as

$$
\nabla W_d(\vec{r})=\frac{\alpha_d}{h^{d+1}}\frac{\vec{r}}{|\vec{r}|}\phi'\left(\frac{|\vec{r}|}{h}\right).$$

Poly6, Spiky, Cubic, Quintic, Gaussian kernels are constructed in this way and are verified.

In `src/common/Kernels.h`, `UnitKernel::Weight()` returns $\alpha_d\phi(|\vec{r}|/h)$, and `UnitKernel::Grad()` returns $\alpha_d\phi'(|\vec{r}|/h)$. And `KernelSPH` class takes care for $h$ and $\vec{r}/|\vec{r}|$.

# Poly6 Kernel
$$
\alpha_d=
\begin{cases}
\frac{35}{32},& d=1\\
\frac{4}{\pi},& d=2\\
\frac{315}{64\pi}, & d=3
\end{cases}
$$

$$
\phi(r)=
\begin{cases}
{(1 - r^2)}^3,& 0 \leq r <1\\
0, & otherwise
\end{cases}
$$

$$
\phi'(r)=
\begin{cases}
{6r(r^2 - 1)^2},& 0 \leq r <1\\
0, & otherwise
\end{cases}
$$

See: Particle-based fluid simulation for interactive applications(2003)


# Spiky Kernel

$$
\alpha_d=
\begin{cases}
2,& d=1\\
\frac{10}{\pi},& d=2\\
\frac{15}{\pi}, & d=3
\end{cases}
$$

$$
\phi(r)=
\begin{cases}
{(1 - r)}^3,& 0 \leq r <1\\
0, & otherwise
\end{cases}
$$

$$
\phi'(r)=
\begin{cases}
-3{(1 - r)}^2,& 0 \leq r <1\\
0, & otherwise
\end{cases}
$$

See: Particle-based fluid simulation for interactive applications(2003)

# Cubic Kernel

Derived from B-Spline.

$$
\alpha_d=
\begin{cases}
\frac{4}{3},& d=1\\
\frac{40}{7\pi},& d=2\\
\frac{8}{\pi}, & d=3
\end{cases}
$$

$$
\phi(r)=
\begin{cases}
6r^3-6r^2+1,& 0 \leq r < 0.5\\
2(1-r)^3, & 0.5 \leq r < 1\\
0, & otherwise
\end{cases}
$$

$$
\phi'(r)=
\begin{cases}
6r(3r-2),& 0 \leq r < 0.5\\
-6(1-r)^2, & 0.5 \leq r < 1\\
0, & otherwise
\end{cases}
$$

See: A refined particle method for astrophysical problems (1985)

# Quintic Kernel

$$
\alpha_d=
\begin{cases}
\frac{1}{40},& d=1\\
\frac{63}{478\pi},& d=2\\
\frac{81}{359\pi}, & d=3
\end{cases}
$$

$$
\phi(r)=
\begin{cases}
(3-3r)^5-6(2-3r)^5+15(1-3r)^5,& 0 \leq r < \frac{1}{3}\\
(3-3r)^5-6(2-3r)^5, & \frac{1}{3} \leq r < \frac{2}{3}\\
(3-3r)^5, &\frac{2}{3}\leq r<1\\
0, & otherwise
\end{cases}
$$

$$
\phi'(r)=
\begin{cases}
-15(3-3r)^4+90(2-3r)^4-225(1-3r)^4,& 0 \leq r < \frac{1}{3}\\
-15(3-3r)^4+90(2-3r)^4, & \frac{1}{3} \leq r < \frac{2}{3}\\
-15(3-3r)^4, &\frac{2}{3}\leq r<1\\
0, & otherwise
\end{cases}
$$

See: Analysis of smoothed particle hydrodynamics with applications (1996)

# Gaussian Kernel

There is an additional parameter $k$ controls the sharpness of the kernel: $\sqrt{2}\sigma=1/k$. In code it's `real UnitGAUSSIAN::trunc_num`, default as $3$.

Gaussian kernel is not really truncated at $1$, but when $k$ is large enough, the function value beyond $1$ will be very close to $0$.

$$
\alpha_d=
\begin{cases}
\frac{k}{\sqrt{\pi}}=\frac{1}{\sigma\sqrt{2\pi}},& d=1\\
\frac{k^2}{\pi^2}=\frac{1}{2\sigma^2\pi},& d=2\\
\frac{k^3}{\pi^3}=\frac{1}{(\sigma\sqrt{2\pi})^3}, & d=3
\end{cases}
$$

$$
\phi(r)=e^{-(kr)^2}=e^{-\frac{r^2}{2\sigma^2}}.
$$

$$
\phi'(r)=-2k^2re^{-(kr)^2}=-\frac{r}{\sigma^2}e^{-\frac{r^2}{2\sigma^2}}.
$$

See: Smoothed particle hydrodynamics: theory and application to non-spherical stars (1977)

# <font color='red'> Deprecated Kernels

Legacy below. Removed in current code.

</font>

# Viscosity Kernel
$$
W_{viscosity}(r,h)=\alpha_d \times  
\begin{cases}
-\frac{r^3}{2h^3}+\frac{r^2}{h^2}+\frac{h}{2r}-1,& 0 \leq r \leq h\\
0, & otherwise
\end{cases}
$$

$$
\nabla_{viscosity}(r,h)=\alpha_d \times \frac{\vec{r}}{|\vec{r}|}
\begin{cases}
-\frac{3r^2}{2h^3}+\frac{2r}{h^2}-\frac{h}{2r^2}-1,& 0 \leq r \leq h\\
0, & otherwise
\end{cases}
$$

$$
    \nabla^2_{viscosity}(r,h)= \frac{45}{\pi h^6}(h-r), d=3
$$

$$
\alpha_d=
\begin{cases}
0,& d=1\\
\frac{10}{3\pi h^2},& d=2\\
\frac{15}{2\pi h^3}, & d=3
\end{cases}
$$

See: Particle-based fluid simulation for interactive applications(2003)

<font color='red'> TODO:What's that? </font>
This is a kernel designed for computation of viscosity force. Currently not implemented in codes.

# Fourth Interpolation Kernel

$$
W_{intp4}(r,h)=\alpha_d
\begin{cases}
{(3-\frac{3r}{h})}^5-6{(2-\frac{3r}{h})}^5+15{(1-\frac{3r}{h})}^5,&0 \leq r < \frac{h}{3}\\
{(3-\frac{3r}{h})}^5-6{(2-\frac{3r}{h})}^5,&\frac{h}{3} \leq r < \frac{2h}{3}\\
{(3-\frac{3r}{h})}^5,&\frac{2h}{3} \leq r < h\\
0, & otherwise
\end{cases}
$$

$$
\nabla_{intp4}(r,h)=\alpha_d \frac{\vec{r}}{|\vec{r}|}
\begin{cases}
-\frac{15}{h}{(3-\frac{r}{h})}^4+\frac{90}{h}{(2-\frac{r}{h})}^4-\frac{225}{h}{(1-\frac{r}{h})}^4,&0 \leq r < \frac{h}{3}\\
-\frac{15}{h}{(3-\frac{r}{h})}^4+\frac{90}{h}{(2-\frac{r}{h})}^4,&\frac{h}{3} \leq r < \frac{2h}{3}\\
-\frac{15}{h}{(3-\frac{r}{h})}^4,&\frac{2h}{3} \leq r < h\\
0, & otherwise
\end{cases}
$$

$$
\alpha_d=
\begin{cases}
,& d=1\\
\frac{63}{478\pi h^2},& d=2\\
\frac{81}{359\pi h^3}, & d=3
\end{cases}
$$

See: A smoothed particle hydrodynamics model for miscible flow in three-dimensional fractures and the two-dimensional Rayleigh--Taylor instability (2005)

# Kernel Laplacians
I believe they're wrong because I derived different Laplacian formulas with MATLAB, and they're numerically errorsome. SPH algorithm does not need Laplacians of kernels, also. However I decided to keep them here, waiting for some mathematician to figure it out.

$$
\nabla^2_{poly6}(r,h)=\alpha_d \times
\begin{cases}
{(r^2 - h^2)(3h^2-7r^2)},& 0 \leq r \leq h\\
0, & otherwise
\end{cases}
$$

$$
\nabla^2_{spiky}(r,h)=\alpha_d \times
\begin{cases}
\frac{-6}{r}(h-r)(h-2r),& 0 \leq r \leq h\\
0, & otherwise
\end{cases}
$$

# Author(s)
Xiangxin Kong, Mengdi Wang
