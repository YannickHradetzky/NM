Let $\Delta$ be a partition of $[a,b]$ consisting of $n+1$ and $y_{0}, \dots, y_{n} \in R$. We want to find a linear interpolating Spline, that is, a function $s \in S^{1}(\Delta)$ satisfying:

$$
s(x_{i})=y_{i},     \text{    for } i=0,\dots,n
$$
The Space $S^{1}(\Delta)$ has dimension $n+1$. An the interpolation problem above imposes $n+1$ conditions. Therefore the linear spline interpolation problem has a unique solution. 

**Using Lagragian polynomials**
Looking at the subinterval $[x_{i},x_{i+1}]$:
$$
l_{i}(x) = y_{i} \frac{x-x_{i+1}}{x_{i}-x_{i+1}}+y_{i+1}\frac{x-x_{i}}{x_{i+1}-x_{i}}
$$

**Using Delta Functions**
Another function satisfying the interpolation conditions are:
$$
\Lambda_{i}(x_{j}) = \delta _{ij} \text{   for }j=0,\dots,n
$$

# Error estimation

## Theorem 3.5 
Let $f \in C^{2}[a,b]$ and denote $s \in S^{1}(\Delta)$ the linear spline that interpolates $f$ at some partition $\Delta$ if $[a,b]$

Applying **Theorem 3.3** in [[Polynomial Interpolation]] to a subinterval $[x_{i-1},x_{i}]$ gives the estimate
$$
\lvert p(x)-f(x) \rvert \leq \frac{\lvert \lvert f^{(n+1)} \rvert  \rvert_{\infty} }{(n+1)!}\prod_{k=0}^{n}\lvert x-x_{k} \rvert 
$$
So in this case with $h_{i}=x_{i}-x_{i-1}$ and $h_{max}=\max_{i}h_{i}$:

$$
\lvert s(x)-f(x) \rvert \leq \frac{\lvert \lvert f'' \rvert  \rvert _{\infty}}{8} h_{max}^{2}
$$
