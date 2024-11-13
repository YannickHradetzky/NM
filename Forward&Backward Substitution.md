#NM 
[[Gaussian Elimination]]

linear equation $Ax=b$, when $A \in K^{n\times n}$ is regular. 
- Direkt Methods: An algorithm  that produces a solution within a finite number of steps
- Iterative Methods: Take infinite number of steps  and have to terminate at some point 

# Direkt Methods for the problem $Ax=b$

## Matrix factorization or Matrix decomposition

- $A=BC$ with $B,C \in K^{n\times n}$ then $B\underbrace{ Cx }_{ y } = b$

### Solving
1) find the Matrizes $B,C$
2) Solve $By=b$ for $y$
3) Solve $Cx=y$ for $x$

### Pitfalls 
A relative error in $b$ will be amplified by the a factor $\kappa(B)$ to produce a new errr in $y$. Similarly, the relative error in $y$ can be magnified by $\kappa(C)$. In total a relative error in $b$ propagates to the final result $x$ with an amplifying factor of $\kappa(B)\kappa(C)$. Therefore solving $Ax=b$ via factorization only makes sense if $\kappa(B)\kappa(C)\leq \kappa(A)$

# Triangular Systems 

## Forward substitution 

Let $L=(l_{ij}) \in K^{n\times n}$ be a lower triangular matrix where $l_{ij}=0$ for $i<j$. 
There is an obvious procedure for solving such a system in $n$ steps. 

### Solving 
1) Solve the first equation $l_{11}x_{1}=b_{1} \Leftrightarrow x_{1}=\frac{b_{1}}{l_{11}}$
2) Plug this in the second equation to solve for $x_{2}$
3) Keep going which leads to the general form 

$$
\begin{align}
x_{k}=\frac{\left( b_{k}-\sum_{j=1}^{k-1}l_{kj}x_{j} \right)}{l_{kk}},    k=2,\dots,n
\end{align}
$$

## Theorems 
### 2.1 
If $L \in K^{n\times n}$ is a regular lower triangle matrix, then for every $b \in K^{n}$ forward substitution yields a unique solution of $Lx=b$
#### Proof
All upper and lower triangular matrixes have non-zero determinants. And therefore no zeros on the diagonal. Therefore it is invertible and forward substitution can not fail. 

### 2.2 
For $L\in K^{n \times n}$ and $b\in K^{n}$ forward substitution requires $n^{2}$ flops

#### Proof 
Think of how many operations are needed in the $k$ step and then sum this over the $n$ equations that need to be solved

### 2.3 
Forward substitution is backward stable. That is, for all lower triangular matrices $L \in F^{n\times n}$ and $b\in F^{n}$ 

## Backward substitution 

For an upper triangular matrix $U \in K^{n \times n}$ where $u_{ij}=0$ if $i<j$ the system $Ux=b$. 
The Algorithm for solving such a system is analogous to the forwards substitution. 

### Solving
1) Solve the last equation $u_{nn}x_{n}=b_{n} \Leftrightarrow x_{n}=\frac{b_{n}}{u_{nn}}$
2) Plug this into the second to last line
3) Repeat 

$$
\begin{align}
x_{k} = \frac{\left( b_{k}-\sum_{j=k+1}^{n}u_{k,j}x_{j} \right)}{u_{k,k}}
\end{align}
$$


