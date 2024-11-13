# 1) Determinant of matrix factorizations

## a) LU factorization

So for a lower triangular matrix:

$$
 L = \begin{pmatrix}
1 & 0 & \dots & 0 \\
l_{21} & 1 & \dots & 0 \\
\vdots & \dots & 1 & 0 \\
l_{n,1} & \dots & \dots & 1
\end{pmatrix}
$$
Therefore $\det L = \prod_{i=1}^{n}l_{ii}$

And for an upper triangle matrix: 

$$
U = \begin{pmatrix}
1 & u_{12} & \dots &  u_{1n} \\
0 & 1 & \dots & u_{2n} \\
\vdots  & \dots & \ddots  & u_{kn} \\
0 & \dots & \dots & 1
\end{pmatrix}
$$
Also here $\det U = \prod_{i=1}^{n}u_{ii}$

So in total $\det A=\det L\cdot \det U=\prod_{i=1}^{n}U_{ii}$

## b) LU factorization with partial pivoting

$$
\begin{align}
PA & =LU \\
\det(PA) & =\det(LU) \\
\overbrace{ \det P }^{ \begin{cases}
1 & \text{if even swaps} \\
-1 & \text{if odd swaps}
\end{cases} }\cdot \det A & = \overbrace{ \det L }^{ 1 }\cdot \det U \\
\det A &  = \frac{\det U}{\det P} = \begin{cases}
1\cdot \det U & \text{if even swaps} \\
-1\cdot \det U & \text{if odd swaps}
\end{cases}
\end{align}
$$
And since $U$ is an upper triangle matrix the determinant is given by its diagonal entries. 

## c) 

Since the matrix needs to be SPD (Self-adjoint positive definite) 

$$
\begin{align}
A  & = R^{*}R \\
\det A& =\det(R^{*}R) \\
& =\det R^{*} \det R \\
& =(\det R)^{*}\det R \\
& = \lvert \det R \rvert ^{2}
\end{align}
$$
## d) 

$$
\begin{align}
\det A  & =\det Q\cdot \det R \\
 & = \det Q \cdot \prod_{i=1}^{n}R_{ii} \\
 & = \lvert \det Q \rvert \cdot \prod_{i=1}^{n}R_{ii} \\
 & =\begin{cases}
1 \cdot\prod_{i=1}^{n}R_{ii} \\
-1\cdot \prod_{i=1}^{n}R_{ii}
\end{cases}
\end{align}
$$
We definitely know that the magnitude of $\det Q$ is $1$. 
$\det Q$ can be a complex number with magnitude $1$. So we might have a certain phase but the magnitude is definitely $1$.



# 3) Show that complexity of LU factorization of band matrix is linear in $n$

## Operations needed for processing row $k$

**Dividing by the Pivot $a_{k,k}$**

Since each row contains a maximum of $2p+1$ non zero Elements, the number of divisions has an upper bound of $2p+1$ 
$$
\text{Number of Divisions in Step } k = 2p+1 
$$
**Eliminate Entries below the diagonal**

Again using the maximum number of non zero entries $2p+1$:
Each row contributes a maximum of $2p+1$ subtractions and $2p+1$ multiplications
Because of the band structure we can put a maximum on the number of rows that need to be eliminated. At each step only rows ranging from $k+1$ to $k+p$ are involved. Because in the other rows the entries are already $0$. This means that there are at most $p$ rows involved.

$$
\text{Number of Steps in Elmination } = p \cdot (4p+2) = 4p^{2}+ 2p
$$
**Total Operations in Step k**

$$
\text{Total Operations in Step }k=2p+1+4p^{2}+2p = 4p^{2}+4p+1
$$

## Total Operations

Since we are processing $n$ rows the total number of flops is bound by: 

$$
\text{Total Amount of Operations} \leq n\cdot (4p^{2}+4p+1)
$$

Which scales linear with $n$

Task 




# 4) Cholesky factorization of: 


$$
A = \begin{pmatrix}
2 & 1 & 0 & 0  \\
1 & 2 & 1 & 0 \\
0 & 1 & 2 & 1 \\
0 & 0 & 1 & 2
\end{pmatrix}
$$

## General Formulas
Given: $$ l_{k,i} = \frac{a_{k,i} - \sum_{j=1}^{i-1} l_{ij} \cdot l_{kj}}{l_{ii}} $$ and $$ l_{kk} = \sqrt{a_{kk} - \sum_{j=1}^{k-1} l_{kj}^2} $$**Step-by-Step Calculation** 
1. Calculate $l_{11}$ $$ l_{11} = \sqrt{a_{11}} = \sqrt{2} $$2. Calculate $l_{21}$ $$ l_{21} = \frac{a_{21}}{l_{11}} = \frac{1}{\sqrt{2}} $$3. Calculate $l_{22}$ $$ l_{22} = \sqrt{a_{22} - l_{21}^2} = \sqrt{2 - \frac{1}{2}} = \sqrt{1.5} = \sqrt{\frac{3}{2}} = \frac{\sqrt{6}}{2} $$4. Calculate $l_{31}$ $$ l_{31} = \frac{a_{31}}{l_{11}} = 0 $$5. Calculate $l_{32}$ $$ l_{32} = \frac{a_{32} - l_{31} \cdot l_{21}}{l_{22}} = \frac{1 - 0}{\frac{\sqrt{6}}{2}} = \frac{2}{\sqrt{6}} = \frac{\sqrt{6}}{3} $$6. Calculate $l_{33}$ $$ l_{33} = \sqrt{a_{33} - l_{32}^2} = \sqrt{2 - \left( \frac{\sqrt{6}}{3} \right)^2} $$ $$ = \sqrt{2 - \frac{6}{9}} = \sqrt{2 - \frac{2}{3}} = \sqrt{\frac{4}{3}} = \frac{2 \sqrt{3}}{3} $$ 7. Calculate $l_{41}$ $$ l_{41} = \frac{a_{41}}{l_{11}} = 0 $$ 8. Calculate $l_{42}$ $$ l_{42} = \frac{a_{42} - l_{41} \cdot l_{21}}{l_{22}} = \frac{0}{\frac{\sqrt{6}}{2}} = 0 $$ 9. Calculate $l_{43}$ $$ l_{43} = \frac{a_{43} - l_{41} \cdot l_{31} - l_{42} \cdot l_{32}}{l_{33}} = \frac{1 - 0 - 0}{\frac{2 \sqrt{3}}{3}} = \frac{3}{2 \sqrt{3}} = \frac{\sqrt{3}}{2} $$ 10. Calculate $l_{44}$ $$ l_{44} = \sqrt{a_{44} - l_{43}^2} = \sqrt{2 - \left( \frac{\sqrt{3}}{2} \right)^2} $$ $$ = \sqrt{2 - \frac{3}{4}} = \sqrt{\frac{5}{4}} = \frac{\sqrt{5}}{2} $$ **Final Result for $L$ The Cholesky factor $L$ is:** $$ L = \begin{pmatrix} \sqrt{2} & 0 & 0 & 0 \\ \frac{1}{\sqrt{2}} & \frac{\sqrt{6}}{2} & 0 & 0 \\ 0 & \frac{\sqrt{6}}{3} & \frac{2 \sqrt{3}}{3} & 0 \\ 0 & 0 & \frac{\sqrt{3}}{2} & \frac{\sqrt{5}}{2} \end{pmatrix} $$
# 5) If $A\in K^{n\times n}$ has a Cholesky factorization then $A$ is SPD

**Symmetric/Hermitian**
So A has a decomposition like $A=R^{*}R$
$$
A^{*} = (R^{*}R)^{*} = R^{*}R = A
$$

**Positive Definiteness**

A Matrix is positive definite if, for any non-zero vector $x \in K^{n}$, we have:
$$
\vec{x}^{*}\vec{A}x>0
$$
$$
\begin{align}
\vec{x}^{*}\vec{A}x  & = \vec{x}^{*}(R^{*}R)\vec{x} \\
  & = (\vec{x}^{*}R^{*})\underbrace{ (\vec{R}x) }_{ \vec{y} } \\
 & =\vec{y}^{*}\vec{y} = \lvert \vec{y} \rvert ^{2}
\end{align}
$$

We know that $R$ is invertible. Therefore $\vec{y}\neq 0$ if $\vec{x}\neq 0$.
And we know for sure that as long as $\vec{x}\neq 0$  that $\vec{y}^{*}\vec{y}>0$


# 6) 

## b) 
So if **$A$ is self-adjoint (symmetric) but not Positive Definite** then the algorithm will break because some square roots might be zero or undefined. Zero is also not allowed because the matrix $R$ needs to have positive diagonal entires. 

If $A$ **is not self-adjoint but positve definite** the factorization will not work because it counts on $A$ having mirrored entries. 

For a simple $2\times 2$ matrix

$$
\begin{pmatrix}
1 & 2 \\
2 & -3 
\end{pmatrix} = \begin{pmatrix}
r_{11}^{2} & r_{11}r_{12} \\
r_{11}r_{12} & r_{12}^{2}+r_{22}^{2}
\end{pmatrix}
$$

It breaks down because you cant solve the 
$$
r_{11}r_{12}
$$
line for two different solutions 


