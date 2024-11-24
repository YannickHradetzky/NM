A polynomial of degree $n$ or less is entirely determined by its $n+1$ coefficients $a_{0}, \dots, a_{n}$ . 

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
