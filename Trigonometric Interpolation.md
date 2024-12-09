A **real** trigonometric polynomial is a function $p: R\to R$ of the form: 
$$
p(x) = a_{0}+\sum_{k=1}^{N}a_{k}\cos kx + \sum_{k=1}^{N}b_{k}\sin kx
$$
---
**Reduction to polynomial Interpolation**
Define $2N+1$ complex coefficients according to
$$
\begin{align}
\gamma_{0} & = a_{0} \\
\gamma_{k} & = \frac{a_{k}-ib_{k}}{2}, \text{ for } 1\leq k\leq N \\
\gamma_{-k} & = \frac{a_{k}+ib_{k}}{2}, \text{ for } 1\leq k\leq N
\end{align}
$$
Then we can use Euler's formula to rewrite the above polynomial: 
$$
p(x) = \sum_{k=-N}^{N}\gamma_{k}e^{ikx}
$$
with $z =e^{ix}$
$$
p(x) = \sum_{k=-N}^{N}\gamma_{k}z^{k} = z^{_{-N}}\sum_{k=0}^{2N}\gamma_{k-N}z^{k}
$$

And the interpolation condition $p(x_{j})= y_{j}$:
$$
\sum_{k=0}^{2N}\gamma_{k_{N}}z_{j}^{k}=z_{j}^{N}y_{j}
$$
---
**Equidistant Interpolation Points**
Assuming we have uniform partitions of the interval $[0,2\pi]$ with 
$$
x_{j}=\frac{2\pi j}{n}, \text{ j=0,...,n-1}
$$

using $w_{n}=e^{2\pi i/n}$
$$
q(x_{k})=\sum_{j=0}^{n-1}c_{j}e^{ijx_{k}}=\sum_{j=0}^{n-1}c_{j}w^{jk}_{n}
$$
Which can be written in matrix-vector notation 
$$
F_{n}\vec{c}=\vec{y}
$$
**Fourier Matrix properties**
- Fourier matrix is **symmetric** ($F_{n}^{T}=F_{n}$) and $n^{-1/2}F_{n}$ is **unitary**
- The solution to the matrix-vector notation can be written as 
	- $c=1/n \cdot \bar{F_{n}}y=F_{n}^{-1}y$

---
**Obtaining Real polynomial**
- If $n$ is odd:
	- $N=\frac{n-1}{2}$
	- $a_{0}=c_{0}$
	- $a_{k}=c_{k}+c_{n-k}$
	- $b_{k}=i(c_{k}-c_{n-k})$
- If $n$ is even
	- $N=\frac{n}{2}$
	- $a_{0}=c_{0}$
	- $a_{k}=c_{k}+c_{n-k}$
	- $b_{k}=i(c_{k}-c_{n-k})$
	- $a_{N}=c_{N}$
	- $b_{N}=0$
