[[Skript Numerics 1.pdf]]

The aim of numerical integration is to approximate the numerical value of definite integrals. 
$$
I[f]=\int_{a}^{b}f(x)dx
$$
In order to numerically integrate a given function $f$ we consider the **quadrature rules** $Q[f]\approx I[f]$

$$
Q[f]=\sum_{i=0}^{n}w_{i}f(x_{i})
$$
- **In General** we say that a quadrature rile has **degree** $k\in N$ if all $p \in P_{k}$ are **integrated exactly**.
- The degree of 
$$
Q[f]=\sum_{i=0}^{n}w_{i}f(x_{i})
$$
	**cannot** be higher than $2n+1$

---
# Midpoint Rule
The simplest example of a quadrature rule: 
$$
Q[f]=(b-a)f\left( \frac{a+b}{2} \right)
$$
It is **exact** for polynomials  $\in P_{1}$

---
# Trapeziodal rule
$$
Q[f]=\frac{b-a}{2}f(a)+\frac{b-a}{2}f(b)
$$


---
# Newton-Cotes Formulas
Are a special class of quadrature rules- Their basic idea is to replace the integrand by the polynomial $p$ that interpolates $f$ at the nodal points $x_{i}$. So $Q[f]:=I[p]$. The polynomial $p$ can be **integrated exactly** by taking the right **linear combinations** of its values at the nodal points.

$$
\begin{align}
I[f]\approx I[p] & =\int_{a}^{b}p(x)dx=\int_{a}^{b}\overbrace{ \left( \sum_{i=0}^{n}f(x_{i}L_{i}(x)) \right)U }^{ \text{Lagrange Interpolation} }dx \\
 & =\sum_{i=0}^{n}f(x_{i})\underbrace{ \int_{a}^{b}L_{i}(x)dx }_{ w_{i} }=:Q[f]
\end{align}
$$

- **The degree** of a Newton-Cotes formula with $n+1$ nodal points, as defined above, is **at least** $n$.
- Is called **closed** if $a$ and $b$ are nodal points. Otherwise its called **open**

- A **closed** Newton-Cotes formula with $n+1$ nodal points has degree:
$$
r = \begin{cases}
n,  & \text{if n is odd} \\
n+1,  & \text{if n is even}
\end{cases}
$$

---
## Simpson Rule
The closed Newton-Cotes formula for $n=2$ is the **Simpson Rule**
$$
Q[f]=\frac{b-a}{6}\left[ f(a)+4f\left( \frac{a+b}{2} \right)+f(b) \right]
$$
- According the general rule for the degree all polynomials in $P_{3}$ are integrated exactly. 

---
## Error estimates
The quadrature error: 
$$
\lvert Q[f]-I[f] \rvert 
$$
An important tool for obtaining error estimates is the **triangle inequality for integrals**
$$
\left\lvert  \int_{a}^{b}g(x)dx  \right\rvert \leq \int_{a}^{b}\lvert g(x) \rvert dx
$$

**In general**:
A Newton-Cotes formula $Q$ of degree $r$ satisfies
$$
\lvert Q[f]-I[f] \rvert\leq C \frac{\lvert \lvert f^{(n+1)} \rvert  \rvert_{\infty}}{(r+1)!} (b-a)^{r+2}
$$
The value of the constant $C>0$ depends on the distribution of the nodal points. 
It suffers from the **same problems as polynomial interpolation**:
A **large number of equidistant nodal points should be avoided**, unless the derivatives of $f$ can be **guaranteed to have a small $\infty$-norms**.  Fortunately, the very idea that overcomes this problem in the case of interpolation also works for numerical integration. 

---
# Composite Rules
Instead of using Newton-Cotes formulas with a large number of nodal points, it is often better to **subdivide** the domain of integration $[a,b]$ and use quadrature rules with fewer points on every subinterval. 

$$
I[f]=\int_{a}^{b}f(x)dx=\sum_{j=1}^{N}\int_{a_{j}}^{b_{j}}f(x)dx \approx \sum_{j=1}^{N}Q^{(j)}[f]=:Q_{N}[f]
$$
The error can be estimated the same way as for non composite Newton-Cotes formula divided by $N^{r+1}$.

---
# # Gaussian Quadrature

Gaussian quadrature is a numerical integration method that approximates the integral of a function $f(x)$ over an interval $[a, b]$ using a weighted sum of function evaluations at specific points, called **nodes**. 

## Key Formula
The integral 
$$
\int_a^b f(x) \, dx
$$
is approximated as:
$$
\int_a^b f(x) \, dx \approx \sum_{i=1}^n w_i f(x_i),
$$
where $x_i$ are the nodes and $w_i$ are the weights.

## Advantages
- **Accuracy:** It is exact for polynomials of degree $2n-1$ or less when $n$ nodes are used.
- **Efficiency:** Requires fewer function evaluations compared to other methods for smooth functions.

## Legendre Polynomials
For integration over $[-1, 1]$, the nodes $x_i$ are the roots of the $n$th Legendre polynomial $P_n(x)$, and weights $w_i$ are computed accordingly:
$$
w_i = \int_{-1}^1 \prod_{j \neq i} \frac{x - x_j}{x_i - x_j} \, dx.
$$

## Applications
- Widely used in numerical analysis, physics, and engineering for solving integrals with high precision.

### Notes:
- For arbitrary intervals $[a, b]$, transform to $[-1, 1]$ using $x = \frac{b-a}{2}\xi + \frac{b+a}{2}$, where $\xi$ is the new variable.
  