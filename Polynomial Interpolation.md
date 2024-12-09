A polynomial of degree $n$ or less is entirely determined by its $n+1$ coefficients $a_{0}, \dots, a_{n}$ . 
# General

**The general problem**

Assume to be given $n+1$ real and distinct data points
$$
(x_{0},y_{0}), \dots, (x_{n},y_{n}) \in R^{2}
$$
We look for a polynomial $p$ of degree $n$ such that: 
$$
p(x_{i})=y_{i} \forall i \in \{ 0,\dots,n \}
$$
which leads to a system of $n+1$ linear equations. 

The matrix representing this system is called **Vandermonde matrix**. The polynomial interpolation problem has a unique solution, if the Vandermonde matrix is regular.

**Choice of Basis**

The choice of basis strongly influences the condition number of the resulting Vandermonde matrix.

- monomials 
- newton polynomials
- Lagrange polynomials

## Error Estimates

Assuming that $y_{i}$ are the values of some function $f$, so that: 
$$
y_{i} = f(x_{i}) 
$$
The question : How well does the interpolating polynomial $p$ approximate $f$ at points $x\neq x_{i}$?

Let $C^{k}[a,b]$ be the space of all functions $f : [a,b] \to R$ which have a continuous $k$-th derivative $f^{(k)}$. 

---
### Theorem 3.3

Let$n \in N$ and $a \leq x_{0}<\dots>x_{n}\leq b$. Suppose that $f\in C^{n+1}[a,b]$ and that $p \in P_{n}$ is the unique interpolating polynomial. 

Then for every $x \in[a,b]$, there is a $\xi \in[a,b]$ such that: 
$$
p(x)-f(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!}\prod_{k=0}^{n}(x-x_{k})
$$
This provides us with a formula for the interpolation error at every point $x \in [a,b]$ but we don't know the value of $\xi$. 

We can obtain a pointwise upper bound for the interpolation error:

$$
\lvert p(x)-f(x) \rvert \leq \frac{\lvert \lvert f^{(n+1)} \rvert  \rvert_{\infty} }{(n+1)!}\prod_{k=0}^{n}\lvert x-x_{k} \rvert 
$$
and a normwise upper bound: 

$$
\rvert\lvert p(x)-f(x) \rvert\rvert_{\infty} \leq \frac{\lvert \lvert f^{(n+1)} \rvert  \rvert_{\infty} }{(n+1)!}\max_{x \in[a,b]}\prod_{k=0}^{n}\lvert x-x_{k} \rvert 
$$
This shows that the interpolation error depends on the way the $x_{i}$ are distributed in $[a,b]$.


**Corallary**

If we can put an upper bound on the maximum value of any derivative : $\lvert \lvert f^{(n)} \rvert \rvert_{\infty}\leq M \forall n \in N$. Then $\lim_{ n \to \infty }\lvert \lvert p_{n}-f \rvert \rvert_{\infty}=0$.

Because the factorial grows faster than than the $\prod$ which has an upper bound of $(b-a)^{n+1}$.

**Choosing the points that minimize $\prod_{k}\lvert x-x_{k} \rvert$  are the so-called Chebyshev Nodes**


